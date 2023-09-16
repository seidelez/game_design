using XLSX
using JuMP
using Gurobi
using CSV
using DataFrames
using Pkg
using Plots
using LinearAlgebra
using IterTools

mutable struct Member
    ID
    P_fix
    flex_load
    flex_load_max
    P_max
    lambda_tarif
    BESS
end

mutable struct Variable_Member
    member
    P
    P_plus
    P_minus
    g_plus
    g_minus
end

mutable struct BESS
    SOC_0
    SOC_max
    SOC
    P_max 
end

mutable struct Variable_BESS
    bess
    SoC
    P_max 
end

mutable struct REC
    members
    payout 
    configuration
    power_fix
    power_virtual
    load_fix
    load_virtual
    lambda_pun
    lambda_prem
    alpha
end

mutable struct Var_REC
    REC
    var_members
    P
    P_shared
end

struct coupled_constraint
    A
    b
    n_dimensions
    n_members
end

T=24

function set_SOC(member::Member)
    SOC = member.BESS.SOC_0*ones(T)

    if member.BESS.P_max == 0
        return
    end

    eta = 1
    for i in 1:(T-1)
        SOC[i+1] = SOC[i] + eta*max(0, member.flex_load[i]) + min(0, member.flex_load[i])/eta
    end
    member.BESS.SOC = SOC
end

function empty_BESS()
   return BESS(0, 0, zeros(T), 0)
end

function g(rec::REC)
    res = 0
    Set_total_load(rec)
    Set_total_power(rec)
    P = rec.power_virtual
    L = rec.load_virtual
    if rec.configuration == "charge only"
        L = rec.power_fix
    end
    SE = min.(L , P)
    res = sum(SE)*rec.lambda_prem - sum((P-L).*rec.lambda_pun)
    return res
end

function g_member(rec::REC, member::Member)
    x_delta = deepcopy(member.flex_load)
    x_plus = max.(0, x_delta)
    x_minus = min.(0, x_delta)
    
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)

    res_plus = x_plus.*SE./(L)
    res_minus = (-x_minus).*SE./(P)
    for i = 1:T
        if L[i] == 0
            res_plus[i] = 0
        end
        if P[i] == 0
            res_minus[i] = 0
        end
    end
    
    res = sum(rec.alpha*res_plus + (1-rec.alpha)*res_minus)*rec.lambda_prem - sum(x_delta.*rec.lambda_pun) 

    return res
end

function set_load_fix(rec::REC)
    total_load = 0*rec.members[1].flex_load 
    for member in rec.members
        total_load .+=  max.(0, member.P_fix)
    end
    rec.load_fix = max.(0, total_load)
    return total_load
end

function Set_total_load(rec::REC)
    total_load = 0*rec.load_fix
    for member in rec.members
        total_load .+=  max.(member.flex_load, 0)
    end
    rec.load_virtual = total_load + rec.load_fix
    return total_load
end

function Set_total_power(rec::REC)
    power_res = zeros(T)
    for member in rec.members
        power_res_sum = max.(-member.flex_load, 0)
        if rec.configuration == "charge only" && member.BESS.Pmax == 0
            power_res_sum = 0*power_res_sum
        end
        power_res .+= power_res_sum
    end
    rec.power_virtual = power_res + rec.power_fix
    return power_res
end

function rec_variational_eq(rec::REC)
end

function rebuild_rec(var_rec::Var_REC)
    rec = var_rec.rec
    var_members = var_rec.var_members
    members = []
    for var_member in var_members
        var_member.member.flex_load = value.(var_member.P)
        push!(members, member)
    end
    rec.members = members
    return rec
end

function display_rec(rec::REC)
    Power = Set_total_power(rec) 
    Load = Set_total_load(rec)
    SE = min.(Load, Power)
    println("load last", rec.load_virtual)
    p1 = Plots.plot(1:T, rec.load_virtual, label ="Load")
    println("P_res", rec.power_virtual)
    Plots.plot!(p1, 1:T, rec.power_virtual, label = " Power")
    Plots.plot!(p1, 1:T, SE)
    display(p1)

    p1 = Plots.plot!( 1:T, [member.flex_load  for member in rec.members], label = "flex load" )
    display(p1)
end

function NE_check(rec::REC)
    res = true
    for member in rec.members
        load_saved = deepcopy(member.flex_load)
        model, var_member = optimal_response(rec, member)
        set_silent(model)
        optimize!(model)
        if objective_value(model) > g_member(rec, member) 
            res = false
        end
        member.flex_load = load_saved
    end
    return res
end

function rec_error(rec1::REC, rec2::REC)

    res = 0

    for (i,member) in enumerate(rec1.members)

        v1 = rec1.members[i].flex_load 
        v2 = rec2.members[i].flex_load 
        v = v1-v2

        #println("v1 ", rec1.members[i].flex_load)
        #println("v2 ", rec2.members[i].flex_load)
        #println("norm ", LinearAlgebra.norm(v1-v2))
        res +=  LinearAlgebra.norm(v)
    end

    #res = g(rec2)
    return res, rec_revenue(rec2), V(rec2), NE_check(rec2), is_local_maxima(rec2)
end

function rec_revenue(rec)
    SE = min.(rec.load_virtual, rec.power_virtual)
    revenue = rec.lambda_prem*sum(SE) 
    revenue_SE = rec.lambda_prem*sum(SE) 
    cost = sum(rec.load_virtual.*rec.lambda_pun) 
    gain = sum(rec.power_virtual.*rec.lambda_pun)
    TOTAL = revenue_SE + gain - cost
    return TOTAL
end

function V(rec::REC)
    SE = min.(rec.load_virtual, rec.power_virtual)
    revenue = rec.lambda_prem*sum(SE) 
    revenue_SE = rec.lambda_prem*sum(SE) 
    cost = sum(rec.load_virtual.*rec.lambda_pun) 
    gain = sum(rec.power_virtual.*rec.lambda_pun)

    if rec.payout == "proportional"
        SE_param = 0.5
    elseif rec.payout == "marginal"
        SE_param = 1
    elseif rec.payout == "shared"
        SE_param = 1/N
    else
        SE_param = 1 
    end
    TOTAL = SE_param*revenue_SE + gain - cost
    return TOTAL
end

function equal_rec(rec::REC, rec_2::REC)
    equal = true
    for i = 1:length(rec.members)
        if rec.members[i].flex_load != rec_2.members[i].flex_load
            equal = false
        end
    end
    return equal
end

function rec_desglose(rec)
    Set_total_load(rec)
    Set_total_power(rec)
    SE = min.(rec.load_virtual, rec.power_virtual)
    revenue_SE = rec.lambda_prem*sum(SE) 
    cost = sum(rec.load_virtual.*rec.lambda_pun) 
    gain = sum(rec.power_virtual.*rec.lambda_pun)
    println(" SE ", revenue_SE, " cost ", cost, " gain ", gain, " TOTAL ", revenue_SE + gain - cost)
end

function is_local_maxima(rec::REC)
    eps = 0.001
    res = g(rec)
    for t in 1:T
        for member in rec.members
            save_x = deepcopy(member.flex_load)[t]
            member.flex_load[t] += eps
            if is_feasible(member)
                if g(rec) > res
                    return false
                end
            end

            member.flex_load[t] = save_x
            member.flex_load[t] -= eps
            if is_feasible(member)
                if g(rec) > res
                    return false
                end
            end

            member.flex_load[t] = save_x

        end
    end
    return true
end

function is_feasible(member::Member)
    for t = 1:T
        if member.BESS.P_max != 0
            set_SOC(member)
            if member.flex_load[t] > member.BESS.P_max || member.flex_load[t] < -member.BESS.P_max || 
             member.SOC[t] > member.SOC_max || member.SOC[t] < 0 || 
             member.SOC[0] != member.SOC_0 || member.SOC[T] >= 0.8*member.SOC_0 
                return false
            end
        else
            if member.flex_load[t] > member.BESS.P_max || member.flex_load[t] < 0 || 
                sum(member.P_fix - member.flex_load) > member.flex_load_max 
                return false
            end
        end
    end
    return true
end

function is_REC_feasible(rec::REC)
    for member in rec.members
        if !is_feasible(member)
            return false
        end
    end

    return true
end

function border_direction(rec::REC)
    candidates = []
    for (i, member) in enumerate(rec.members)
        for t = 1:T          
            new_candidate1 = zeros(size(rec.members, 1), T)
            new_candidate2 = zeros(size(rec.members, 1), T)  
            if member.flex_load[t] == 0
                new_candidate1[i, t] = 1
                new_candidate2[i, t] = -1
            end
            push!(candidates, new_candidate1, new_candidate2)
        end
    end
    return candidates
end

function border_direction_alt(rec::REC)
    candidates = []
    new_candidate = zeros(size(rec.members, 1), T)
    eps = 10^(-3)
    model, rec_centralised, costs = optimal_centralised(rec)
    
    @constraint(model, costs >= V(rec) + eps)
    @objective(model, Max, 0)
    set_silent(model)
    optimize!(model)
    println( "model is ", termination_status(model) )

    if termination_status(model) != "OPTIMAL"
        return candidates
    end
    
    for (i, member) in enumerate(rec.members), t in 1:T
        new_candidate[i] = rec_centralised.members[i].flex_load - member.flex_load
    end
    push!(candidates, new_candidate)
    return candidates
end

function border_direction_alt1(rec::REC)
    candidates = []
    new_candidate = zeros(size(rec.members, 1), T)
    eps = 10^(-3)
    model, rec_centralised, costs = optimal_centralised(rec)
    
    @constraint(model, costs >= V(rec) + eps)
    @objective(model, Max, 0)
    set_silent(model)
    optimize!(model)
    println( "model is ", termination_status(model) )

    if termination_status(model) != "OPTIMAL"
        return candidates
    end
    
    for (i, member) in enumerate(rec.members), t in 1:T
        new_candidate[i] = rec_centralised.members[i].flex_load - member.flex_load
    end
    push!(candidates, new_candidate)
    return candidates
end

function border_direction_alt2(rec::REC)
    candidates = []
    for (i, member) in enumerate(rec.members)
        for constraint in vector_constraints(member)
            if constraint[1]'*member.flex_load == constraint[2]
                vectorized_contraint = zeros(size(rec.members, 1), T)
                vectorized_contraint[i] = constraint[1]
                push!(vectorized_contraint, candidates)
            end
        end
    end

    mean_candidate = zeros(size(rec.members, 1), T)
    for candidate in candidates
        mean_candidate += candidate
    end

    mean_candidate /= size(mean_candidate, 1)
    return candidate
end

function vector_constraint(member::Member)
    constraints = []
    for t in 1:T
        c1 = zeros(T)
        c2 = zeros(T)
        c1[t] = 1
        c2[t] = 1
        Pmin = -member.BESS.P_max
        Pmax = max(member.BESS.P_max, member.P_max)
        push!(constraints, [c1,  Pmin])
        push!(constraints, [-c2, Pmax])
    end

    if member.BESS.P_max != 0
        for t in 1:T
            c1 = zeros(T)
            c1[1:t] = ones(t)
            push!(constraints, [c1,  -member.SOC_0])
            push!(constraints, [-c1, -member.SOC_max + member.SOC_0])
        end
    else
        c1 = ones(T)
        push!(constraints, [c1,  sum(member.P_fix)])
    end
    return constraints
end