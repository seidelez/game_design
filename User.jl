using XLSX
using JuMP
#using Gurobi
using CSV
using DataFrames
using Pkg
using Plots
using LinearAlgebra
using IterTools

mutable struct Interval
    a
    b
end

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
    g
    dg
    C
    V
    single_update
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

function interval(a, b)
    ini = min(a, b)
    fin = max(a, b)
    return Interval(a, b)
end


function intersection(I::Interval, J::Interval)
    a = max(I.a, J.a)
    b = min(I.b, J.b)
    return Interval(a, b)
end

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

function summary_rec(rec::REC)
    for (i, member) in enumerate(rec.members)
        println("Member ", i, ": Flex ", member.flex_load_max, ", sum", sum(member.flex_load))
    end
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

function g_member(rec::REC)
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
    total_load = zeros(T)
    for member in rec.members
        total_load .+=  max.(member.flex_load, 0)
    end
    rec.load_virtual = total_load 
    #println("Load virtual ", rec.load_virtual)
    return total_load
end

function Set_total_power(rec::REC)
    power_res = zeros(T)
    for member in rec.members
        power_res = power_res + max.(-member.flex_load, 0)
    end
    rec.power_virtual = power_res + rec.power_fix
    #println("power virtual ", rec.power_virtual)
    return power_res
end

function sigma(rec::REC)
    Set_total_load(rec)
    Set_total_power(rec)
    res = rec.load_virtual .<= rec.power_virtual
    return res
end

function v_power(t1, t2, eta = 1)
    v = zeros(T)
    v[t1] = -1
    v[t2] = +1
    return v
end

function g_cost(rec::REC, member::Member)
    res = -rec.lambda_pun.*member.flex_load
    return sum(res)
end
    
function g_marginal(rec::REC, member::Member)    
    Set_total_load(rec)
    Set_total_power(rec)
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)

    x_plus = max.(0, member.flex_load)
    x_minus = min.(0, member.flex_load)
    SE_without = min.(L - x_plus, P + x_minus)
    
    res =  rec.lambda_prem*(SE - SE_without) 
    return sum(res)
end

function g_marginal_vect(rec::REC, member::Member)    
    Set_total_load(rec)
    Set_total_power(rec)
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)

    x_plus = max.(0, member.flex_load)
    x_minus = min.(0, member.flex_load)
    SE_without = min.(L - x_plus, P + x_minus)
    
    res =  rec.lambda_prem*(SE - SE_without) 
    return res
end

function sum_of_g(gs)
    function g_sum(rec::REC, member::Member)
        res = 0
        for g in gs
            res += g(rec, member)
        end
        return res
    end
    return g_sum
end

function sum_of_dg(dgs)
    function dg_sum(rec::REC, i, j, t1, t2)
        res = 0
        for dg in dgs
            res += dg(rec, i, j, t1, t2)
        end
        return res
    end
    return dg_sum
end

function intersect_of_C( Cs ) 
    function inter_C(rec::REC, member::Member, t1, t2)
        res = Cs[1](rec,member,t1,t2)
        for c in Cs
            res = intersection(res, c(rec,member,t1,t2))
        end
        return res
    end
    return inter_C
end

function dg_cost(rec::REC, i, j, t1, t2)
    res = -rec.lambda_pun[t1] + rec.lambda_pun[t2]
    return res
end

function dg_equal(rec::REC, i, j, t1, t2)
    res = dg_marginal(rec::REC, i, i, t1, t2)
    return res/N
end

function dg_marginal(rec::REC, i, j, t1, t2)
    #IF G NOT DERIVABLE WE TAKE DG_T1^+ AND  DG_T2^- 
    member = rec.members[i]
    member2 = rec.members[j]
    Set_total_load(rec)
    Set_total_power(rec)
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)

    x_plus = max.(0, member.flex_load)
    x_minus = min.(0, member.flex_load)
    SE_without = min.(L - x_plus, P + x_minus)
    SE_marg = (SE - SE_without)

    if i != j
        dg_t1 =  -rec.lambda_prem*((SE_marg .> 0).*(L.>= P))[t1] 
        dg_t2 =  -rec.lambda_prem*((SE_marg .>= 0).*(L.> P))[t2] 
    else
        dg_t1 =  rec.lambda_prem*(L .< P)[t1] 
        dg_t2 =  rec.lambda_prem*(L .<= P)[t2] 
    end
    res = dg_t1 - dg_t2
    return res
end

function dg_penalty(rec::REC, i, j, t1, t2, penalty = 1)
    #IF G NOT DERIVABLE WE TAKE DG_T1^+ AND  DG_T2^- 
    member = rec.members[i]
    member2 = rec.members[j]
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)

    x_plus = max.(0, member.flex_load)
    x_minus = min.(0, member.flex_load)
    SE_without = min.(L - x_plus, P + x_minus)
    SE_marg = (SE - SE_without)
    penalty = x-SE_marg

    if i != j
        dg_t1 =  -rec.lambda_prem*((SE_marg .> 0).*(L.>= P))[t1] - rec.lambda_pun[t1]
        dg_t2 =  -rec.lambda_prem*((SE_marg .>= 0).*(L.> P))[t2] - rec.lambda_pun[t2]
    else
        dg_t1 =  rec.lambda_prem*(L .< P)[t1] - rec.lambda_pun[t1]
        dg_t2 =  rec.lambda_prem*(L .<= P)[t2] - rec.lambda_pun[t2]
    end
    res = dg_t1 - dg_t2
    return res
end

function g_proportional(rec::REC, member::Member)
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)

    x_plus = max.(0, member.flex_load)
    x_minus = min.(0, member.flex_load)

    L_greater = L.>=SE 
    P_greater = P.>=SE 
    
    replace!(L, 0=>1)
    replace!(P, 0=>1)
    res =  rec.lambda_prem*(x_plus.*SE./P + x_minus.*SE./L) 
    return sum(res)
end

function g_equal(rec::REC, member::Member)
    Set_total_power(rec)
    Set_total_load(rec)
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)

    x_plus = max.(0, member.flex_load)
    x_minus = min.(0, member.flex_load)

    res =  rec.lambda_prem*(SE)/N 
    return sum(res)
end

function g_equal_aux(rec::REC, member::Member, g_aux =g_marginal)
    Set_total_power(rec)
    Set_total_load(rec)
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)
    rest = sum(rec.lambda_prem*(SE))
    for member_i in rec.members
        rest -= g_aux(rec, member_i)
    end
    return rest/N
end

function g_equal_aux2(rec::REC, member::Member, g_aux = g_marginal_vect)
    Set_total_power(rec)
    Set_total_load(rec)
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)
    rest_t = rec.lambda_prem*(SE)

    x_plus = max.(0, member.flex_load)
    x_minus = min.(0, member.flex_load)
    SE_without = min.(L - x_plus, P + x_minus)
    N_t = zeros(T)

    for member_i in rec.members
        rest_t -= g_aux(rec, member_i)
        N_t += (member_i.flex_load .> 0)
    end

    rest_t = rest_t.*(member.flex_load.>0)
    replace!(N_t, 0 =>1)
    rest_t = rest_t./N_t
    res = sum(rest_t)
    return res
end

function g_penalty(rec::REC, member::Member, penalty_coef = 1)
    Set_total_load(rec)
    Set_total_power(rec)
    P = rec.power_virtual
    L = rec.load_virtual 
    SE = min.(L , P)

    x_plus = max.(0, member.flex_load)
    x_minus = min.(0, member.flex_load)
    SE_without = min.(L - x_plus, P + x_minus)
    SE_marg = SE - SE_without
    
    
    penalty = rec.lambda_prem*(x_plus .- SE_marg).*(x_plus .>= 0)
    SE_incentive = g_marginal_vect(rec, member) - penalty
    SE_incentive = max.(SE_incentive, 0)
    res = SE_incentive 
    return sum(res)
end

function g_regularized(rec::REC, member::Member, old_member::Member)
    regularizer = LinearAlgebra.norm(member.flex_load - old_member.flex_load, 1)
    return regularizer
end

function dg_proportional(rec::REC, member::Member, t1 , t2)
    x = deepcopy(member.flex_load)
    Power_wo_member = Set_total_load(rec) - max.(0, -x)
    Load_wo_member = Set_total_power(rec) - max.(0, x)
    ordered = t1 < t2
    
    res_cost = (rec.lambda_pun[t1] - rec.lambda_pun[t2])
    function deriv(delta)
        x[t1] -= delta
        x[t2] += delta

        L = Load_wo_member + max.(0, x)
        P = Power_wo_member + max.(0, -x)
        SE = min.(L, P)
        if x[t1] >= 0
            if L <= P
                res += -1
            else
                res += -SE[t1]*Load_wo_member[t1]/(Load_wo_member[t1] + x[t1])^2
            end
        else
            if L >= P
                res += 1
            else
                res += SE[t1]*Power_wo_member[t1]/(Power_wo_member[t1] - x[t1])^2
            end
        end

        if x[t2] >= 0
            if L <= P
                res += 1
            else
                res += SE[t2]*Load_wo_member[t2]/(Load_wo_member[t2] + x[t2])^2
            end
        else
            if L >= P
                res += -1
            else
                res += -SE[t2]*Power_wo_member[t2]/(Power_wo_member[t2] - x[t2])^2
            end
        end

        res = res*sign(delta)
        res = res + res_cost

        return res     
    end 
    
    return deriv
end

function delta_feasible(rec::REC, member::Member, t1, t2, positive = true)
    saved_load = deepcopy(member.flex_load)
    updated_load = deepcopy(member.flex_load)
    member.flex_load = saved_load
    set_SOC(member)
    P_t1 =  updated_load[t1]
    P_t2 =  updated_load[t2]
    P_max = member.P_max 
    P_min = 0
    if member.BESS.P_max != 0
        P_min = -member.BESS.P_max
    end

    Delta_max = min(P_t1 - P_min, P_max - P_t2)
    Delta_min = max(P_t1 - P_max, P_min - P_t2)

    if member.BESS.P_max > 0
        a = minimum(member.BESS.SOC[(t1+1):(t2)])
        b = (maximum(member.BESS.SOC[(t1+1):(t2)])-member.BESS.SOC_max)
        Delta_max = min(Delta_max, a)
        Delta_min = max(Delta_min, b)
        Delta_min = max(Delta_min, 0)
    end

    delta_feasible = Interval(Delta_min, Delta_max)
    return delta_feasible
end

function delta_pos(rec::REC, member::Member, t1, t2)
    return Interval(0, 1e6)
end

function delta_SE(rec::REC, member::Member, t1, t2)
    delta_feas = delta_feasible(rec, member, t1, t2)
    Set_total_power(rec)
    Set_total_load(rec)
    D = rec.load_virtual - rec.power_virtual 
    delta = [D[t2], -D[t1]]
    sort(delta)

    if 0 < delta[1] 
        delta = Interval(0, delta[1])
    elseif  0 < delta[2]
        delta = Interval(delta[1], delta[2])
    else
        delta = Interval(delta[2], 10e5)
    end

    delta = intersection(delta, delta_feas)
    return delta
end

function delta_continious(rec::REC, member::Member, t1, t2, epsilon = 1e-4)
    delta_feas = delta_feasible(rec, member, t1, t2)
    Set_total_power(rec)
    Set_total_load(rec)

    delta_pos = delta_SE(rec, member, t1, t2) 
    delta_neg = delta_SE(rec, member, t2, t1) 

    Delta = Interval(-delta_neg[2], delta_pos[2])
    max_length = 0.33*member.P_max

    if length > max_length
        neg_proportion = abs(Delta.a/length)
        pos_proportion = abs(Delta.b/length)
        Delta = [-neg_proportion*max_length, -pos_proportion*max_length]
    end
    
    if delta_feas.a <0 
        member_neg = deepcopy(member)
        member_neg.flex_load = member_neg.flex_load + v_power(t2, t1)*epsilon
        delta_epsilon_neg = delta_SE(rec, member_neg, t2, t1)
        diff = (delta_neg.a - delta_epsilon_neg.a)^2 + (delta_neg.b - delta_epsilon_neg.b)^2
        if res > eps
            return Delta 
        end
    elseif delta_feas.b >0 
        member_pos = deepcopy(member)
        member_pos.flex_load = member_neg.flex_load + v_power(t1, t2)*epsilon
        delta_epsilon_pos = delta_SE(rec, member_pos, t1, t2)
        diff = (delta_neg.a - delta_epsilon_pos.a)^2 + (delta_neg.b - delta_epsilon_pos.b)^2
        if res > eps
            return Delta 
        end
    elseif !delta.a && !delta.b 
        return Delta
    end

    member_neg = deepcopy(member)
    member_neg.flex_load = member_neg.flex_load + v_power(t1, t2)*Delta.a
    member_pos = deepcopy(member)
    member_pos.flex_load = member_neg.flex_load + v_power(t1, t2)*Delta.b

    return intersection(delta_continuous(rec, member_neg, t1, t2), delta_continious(rec, member_neg, t1, t2) )
end

function delta_postive()
    res = Interval(0, 1e9)
    return res
end

function delta_limited(rec::REC, member, t1, t2)
    eps = 0.001
    power_save = member.flex_load
    member.flex_load = power_save + v_power(t1, t2)*eps
    sigma_plus = sigma(rec)
    member.flex_load = power_save - v_power(t1, t2)*eps
    sigma_minus = sigma(rec)
    member.flex_load = power_save 

    Set_total_load(rec)
    Set_total_power(rec)

    if sigma_plus == sigma_minus
        I1 = Interval(0, sigma_plus[t1]*(rec.power_virtual[t1] - rec.load_virtual[t2]) + (1-sigma_plus[t1])*1e5) 
        I2 = Interval(0, (1-sigma_plus[t2])*(rec.load_virtual[t2] - rec.power_fix[t2]) + (sigma_plus[t1])*1e5) 
        I = intersection(I1, I2)
        res = intersection(I, delta_feasible(rec, member, t1, t2))
    else
        sigma_plus = sigma(member)
        member.flex_load = power_save + v_power(t1, t2)*2*eps
        I = delta_limited(rec, member, t1, t2)
        I.b = I.b + 2*eps
        res = I
        member.flex_load = power_save 
    end

    sum_derivs = 0
    for (j, member2) in enumerate(rec.members)
        sum_derivs += rec.dg(rec, member.ID, j, t1, t2)
    end
    lyap_condition = (sum_derivs >= 0)
    res.b = res.b*lyap_condition 
    return res
end

function V_SE(rec::REC)
    res = 0
    Set_total_load(rec)
    Set_total_power(rec)
    P = rec.power_virtual
    L = rec.load_virtual
    if rec.configuration == "charge only"
        L = rec.power_fix
    end
    SE = min.(L , P)
    res = sum(SE)*rec.lambda_prem 
    return res
end

function V_cost(rec::REC)
    res = 0
    Set_total_load(rec)
    Set_total_power(rec)
    P = rec.power_virtual
    L = rec.load_virtual
    res = -sum((L-P).*rec.lambda_pun)
    return res
end

function V_hybrid(rec::REC, alpha = 0.5)
    return V_cost(rec) + V_SE(rec)
end

function V_g(rec::REC)
    res = 0
    for member in rec.members
        res += rec.g(rec, member)
    end
    return res
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
    verbose = false
    res = true
    eps = 0.001
    for (i, member) in enumerate(rec.members)
        rec2 = deepcopy(rec)
        member2 = rec2.members[i]
        if verbose
            println("flex total ini ", sum(member.flex_load))
            println("Initial g ", rec.g(rec, member), " ")
            println("P ini ", member.flex_load)
        end
        model, objective, var_member = optimal_response(rec2, member2)
        if !verbose
            set_silent(model)
        else
            println("P fin ", member.flex_load)
        end
        optimize!(model)
        objective = objective_value(model)
        member2.flex_load = value.(var_member.P)
        if verbose
            println("P BR ", member2.flex_load)
            println("flex total end ", sum(member.flex_load))
            println("Final ", rec2.g(rec2, member2))
            println("Final obj ", objective)
        end
        if rec.g(rec, member) - rec2.g(rec2, member2) > 0.1
            if verbose            
                println("problemas en ", i)
                println("objective ", objective)
                println("member ", member.flex_load)
                println("member2 ", member2.flex_load)
                println("member2 var", value.(var_member.P))
            end
        end
        if rec2.g(rec2, member2) - rec.g(rec, member)> eps
            return false
        end
    end
    return res
end

function rec_error(rec1::REC, rec2::REC)

    res = 0

    for (i,member) in enumerate(rec1.members)

        v1 = rec1.members[i].flex_load 
        v2 = rec2.members[i].flex_load 
        v = v1-v2
        res +=  LinearAlgebra.norm(v)
    end

    #res = g(rec2)
    return res, rec_revenue(rec2), rec2.V(rec2), NE_check(rec2), is_local_maxima(rec2)
end

function rec_revenue(rec)
    Set_total_load(rec)
    Set_total_power(rec)
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
    println(" SE ", revenue_SE, " cost ", cost, " gain ", gain, " TOTAL ", rec.V(rec))
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