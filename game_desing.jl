using XLSX
using JuMP
#using Gurobi
using CSV
using DataFrames
using Pkg
using Plots, PlotlyJS
using LinearAlgebra
using Roots
using Printf
using JLD2
include("User.jl")
include("bargaining.jl")
include("proximal_gradient_descent.jl")
include("import_data.jl")
include("lyapunov.jl")

Base.show(io::IO, f::Float64) = @printf(io, "%1.3f", f)

function display_rec_load(rec)
    df = DataFrame(timestep=Float64[], ID=Float64[], load=Float64[])
    for member in rec.members
        for t in 1:T
            push!(df, [t, member.ID, member.flex_load[t]])
        end
    end
    p = PlotlyJS.plot(df, x=:timestep, y=:load, color=:ID, kind="bar", Layout(title="Load time series", barmode="relative"))
    display(p)

end

#Payout Function for individual player for strategy s_i
function payout_proportional(member, rec)
    #Water fill to maximize energy 
    member.flex_load = 0*member.flex_load
    s_minus_i = Set_total_load(rec)
    load_initial = s_minus_i
    load_unalocated = member.flex_load_max

    for i in 1:length(member.flex_load)
        excess_res = rec.power_virtual[i] - s_minus_i[i]
        member.flex_load[i] = max(min(excess_res, load_unalocated),0)
        load_unalocated -= member.flex_load[i]
    end
       
    if load_unalocated <= 0
        return 
    end

    k = ( sum(s_minus_i) + member.flex_load_max )/sum(rec.power_virtual)
    member.flex_load = k*rec.power_virtual - s_minus_i
end

function payout_g2(member::Member, rec::REC)
    #candidates dictionnaire of possible candidates
    candidates = []
    member.flex_load = 0*member.flex_load
    s_minus_i = Set_total_load(rec)
    load_unalocated = member.flex_load_max

    for i in 1:length(member.flex_load)
        excess_res = rec.power_virtual[i] - s_minus_i[i]
        member.flex_load[i] = max(min(excess_res, load_unalocated),0)
        load_unalocated -= member.flex_load[i]
    end

    #WE NEED TO UPDATE REC.LOAD_virtual
    fixed_flex_load = member.flex_load
    Set_total_load(rec)

    function payout(x)
        L = rec.power_virtual
        P = rec.load_virtual
        K = max.(P-L, 0*P)
        res = x*(L + min(K , x))/(x+L)
        return res
    end

    function g_alpha_ex(a, except)
        x = zeros(T)
        for t in 1:T
            l = rec.power_virtual[t]
            p = rec.load_virtual[t]
            if l >= p
                x[t] = (-l*sqrt(a) + l)/sqrt(a)
            else
                x[t] = (-l*sqrt(a) + sqrt(l*p))/sqrt(a)
            end
        end

        for t in except
            x[t] = 0
        end
        return sum(x) - load_unalocated     
    end

    function rebuild_x(rec, a, except)
        res = zeros(T)
        for t in 1:T
            l = rec.power_virtual[t]
            p = rec.load_virtual[t]
            if l >= p 
                res[t] = (-l*sqrt(a) + l)/sqrt(a)
            else
                res[t] = (-l*sqrt(a) + sqrt(l*p))/sqrt(a)
            end
        end

        for i in except
            res[i] = 0
        end
        return res
    end

    if load_unalocated <= 0
        return
    end

    print("load_unalocated ", load_unalocated)

    function  g_alpha(x) 
        g_alpha_ex(x, []) 
    end

    alpha_sol = find_zeros(g_alpha, (0.001,  0.99))
    a = alpha_sol
    println("alpha ", a)

    times_power_is0 = 12
    println("member ID ", member.ID, " bo ", a, " ", length(a))
    if length(a) != 0
        member.flex_load = rebuild_x(rec, a[1], [])
    end
  
    #model, var_member = optimal_response(rec, member)
    optimize!(model)
    println( "obj value ", objective_value(model) )
    P_plus = max.(value.(var_member.P), 0)
    P_minus = max.(-value.(var_member.P), 0)
    println( "G_plus", value.(var_member.g_plus))
    println( "G_minus", value.(var_member.g_minus))
    member.flex_load = value.(var_member.P)
    
    println("member ", member.flex_load)
end

function payout_g(member::Member, rec::REC)
    member.flex_load = 0*member.flex_load
    Set_total_load(rec)

    function payout(x)
        L = rec.power_virtual
        P = rec.load_virtual
        K = max.(P-L, 0*P)
        res = x*(L + min(K , x))/(x+L)
        return res
    end

    model, var_member = optimal_response(rec, member)
    set_silent(model)
    optimize!(model)
    results_summary(model, rec, member, var_member)

    return model, var_member
end

function payout_g2(member, rec, g)
    #candidates dictionnaire of possible candidates
    function lagrangian(g, member, rec, T)
        function subdiff_g(x)
            x0 = x[1:T]
            xbar = x[T+1]
            res = zeros(T)
            eps = 1e-3
            for i in 1:T
                e_i = zeros(T)
                e_i[i] =  1
                res[i] = (g(x0+eps*e_i, rec)-g(x0-eps*e_i, rec))/(2*eps)
            end

            res[i] += xbar.*ones(T)
            res += 10*broadcast(abs, member.flex_load_max.*ones(T)-sum(x0))
        end
        return subdiff_g
    end

    #candidates = Dict{Vector::Float64, Float64}()
    subdiff_g = lagrangian(g, member, rec, T)
    lagrangian_sol = Roots.find_zeros(subdiff_g, x0 = [member.flex_load, 1])
    sol = lagrangian_sol[1:T]
    print(sol)
    #candidate[sol] = g(sol)

    for key in candidates.keys()
        if candidates[key] > rec.payout(member.flex_load)
            member.flex_load = key
        end
    end
end

function best_response_i(member, rec, payout = "proportional")

    function func_g(x::Vector{Float64}, rec::REC)
        res = 0
        for t in 1:T
            g_i = x[t]*(rec.power_virtual[t])/(x[t] + rec.load_virtual[t])
            res = res + g_i 
        end
        return res
    end

    if payout == "proportional"
        println("before ",member)
        payout_proportional(member, rec)
        println("after ",member)
        return
    else
        payout_g(member, rec)
    end
end

function best_response_dynamics(rec, bess = "joint", upload = true, learning_rate = 1, verbose = false)
    n_iter = 10
    error = []
    for i in 1:n_iter
        rec_t = deepcopy(rec)
        new_members = []
        for member in rec.members
            saved_load = deepcopy(member.flex_load)
            aux_load = deepcopy(member.flex_load)
            model, var_member = optimal_response(rec, member)
            

            if bess == "separated" && member.BESS.P_max!=0
                model_plus, var_member_plus = optimal_response(rec, member)
                @constraint(model_plus, [i=1:T], var_member_plus.P_minus[i] == min(0, aux_load[i]) )
                set_silent(model_plus)
                optimize!(model_plus)

                for i in 1:T
                    if aux_load[i] >= 0 
                        aux_load[i] = value(var_member_plus.P_plus[i])
                    end
                end

                model_minus, var_member_minus = optimal_response(rec, member) 
                @constraint(model_minus, [i=1:T], var_member_minus.P_plus[i] == max(0, aux_load[i]) )
                set_silent(model_minus)
                optimize!(model_minus)

                for i in 1:T
                    if aux_load[i] <= 0 
                        aux_load[i] = value(var_member_minus.P_minus[i])
                    end
                end
                
                if verbose
                    results_summary(model_minus, rec, member, var_member_minus)
                end
                P_new = aux_load
                model = model_plus
            else
                set_silent(model)
                optimize!(model)
                if verbose
                    esults_summary(model, rec, member, var_member)
                end
                P_new = value.(var_member.P)
            end

            if upload
                member.flex_load = learning_rate*P_new + (1-learning_rate)*saved_load
            else
                new_member = deepcopy(member)
                new_member.flex_load = learning_rate*P_new + (1-learning_rate)*saved_load
                push!(new_members, new_member)
            end
        end

        if !upload
            rec.members = new_members
        end 
        push!(error,rec_error(rec_t, rec))
    end
    return error
end

function upload_rec(rec::REC)
    n_iter = 10
    n_iter_member = 10

    for i in n_iter
        for member in rec.members
            for j in n_iter_member
                model, g = update_member(rec, member)
            end

            println("flex_load ", member.flex_load)
            println("incentive ", g)
        end
    end
end

function results_summary(model, rec, member, var_member)
    lambda_pun = rec.lambda_pun
    if member.BESS.SOC_max != 0
        println("STORAGE")
    end

    if !isa(model, Int)       
        println( "obj value ", objective_value(model) )
        println( "model is ", termination_status(model) )
    end
    println( "Total Load ", rec.load_virtual )
    println( "Total power ", rec.power_virtual )
    P_plus = max.(value.(var_member.P), 0)
    P_minus = max.(-value.(var_member.P), 0)
    println( "costs ", sum(lambda_pun.*var_member.P) )
    println( "G_plus ", value.(var_member.g_plus)*lambda_prem)
    println( "G_minus ", value.(var_member.g_minus)*lambda_prem)
    
    println("member ", value.(var_member.P))
    println("p plus ", value.(var_member.P_plus))
    println("p minus ", value.(var_member.P_minus))
    println("SE incentives ", sum(value.(var_member.g_plus + var_member.g_minus)*lambda_prem))
    println("Costs ", sum(value.(var_member.P).*lambda_pun))
end

members = []
T = 24
Pmax = 5
SOC_max = 2*Pmax
lambda_pun = lambda_pun[1:T]
lambda_prem = 0.03
println("lambda_pun ", lambda_pun)
println("lambda_prem ", lambda_prem)
for i in 1:6
    if i <= 4
        bess = empty_BESS()
        P_fix = 3*ones(T)
        new_member = Member(i, P_fix, P_fix, 0.5*sum(P_fix), Pmax, lambda_pun, bess)
    else
        bess = BESS(0, SOC_max, zeros(T), Pmax)
        new_member = Member(i, zeros(T), zeros(T), Pmax, Pmax, lambda_pun, bess)
    end
    push!(members, new_member)
end

members2 = []
for i in 1:8
    if i in [3,4,5,6]
        bess = empty_BESS()
        P_fix = 0.5*Pmax*ones(T) 
        new_member = Member(i, P_fix, P_fix, sum(P_fix), Pmax, lambda_pun, bess)
    elseif  i in [0] 
        P = 5*ones(T)
        for t = 1:T
            if t%2 == 1 

                P[t] = -P[t]
            end
        end
        bess = BESS(0, 2*SOC_max, P, Pmax)
        new_member = Member(i, zeros(T), 0*P, 90, 5, lambda_pun, bess)
    else 
        bess = empty_BESS()
        P_fix = 5*ones(T)
        new_member = Member(i, P_fix, P_fix, 0, Pmax, lambda_pun, bess)
    end
    push!(members2, new_member)
end
members2 = members2[1:8]

power = 20*rand(Float64, T)
for i in 1:length(power)
    power[i] = 3*min(i, 24-i)
    if (i <= 6) || (i>=18)
        power[i] = 0
    end
end

println("power", power)

power = 2*power
lambda_pun = 100  *lambda_pun

lambda_pun[16:24] = lambda_pun[16:24]
lambda_prem = 0.12
lambda_prem = 1*lambda_prem
payout = "marginal"
configuration = "standalone"
alpha = 0.3
load = ones(T)
single_update = false
Ve = V_hybrid

C = intersect_of_C([delta_feasible, delta_pos])
g2 = sum_of_g([g_cost, g_marginal])
dg = sum_of_dg([dg_cost, dg_marginal])


rec = REC(members2, payout, configuration, power, power, load, load, lambda_pun, lambda_prem, alpha, g2, dg, C, Ve, single_update)
Set_total_load(rec)
Set_total_power(rec)
rec_centralised = deepcopy(rec)
rec_initial = deepcopy(rec)
println("load", rec.load_virtual)
println("P_res", rec.power_virtual)

function reset_rec(rec::REC)
    rec = rec_initial
    return rec
end

#BESTRESPONSEDYNAMICS
type = "proportional"
type = "no"
N = length(rec.members)

rec.payout = "marginal"
print("c", lambda_pun)
err = 0
err, rec_timeseries = lyapunov(rec)
rec = rec_timeseries[end]
println("V(rec)", rec.V(rec))
display_rec(rec)
display_rec_load(rec_timeseries[end])
rec_desglose(rec_timeseries[end])
println(" Error 1", length(err))
println(" Error 1", err)

for member in rec.members
    for t1 in 1:T
        for t2 in 1:T
            if rec.dg(rec, member.ID, member.ID, t1, t2) >= 0 && delta_feasible(rec, member, t1, t2).b > 0
                println("dg ", member.ID, ", t1 ", t1, ", t2 ", t2, " =: ", rec.dg(rec, member.ID, member.ID, t1, t2) )
                println("delta limited " , rec.C(rec, member, t1, t2))
                println("delta feasible ", delta_feasible(rec, member, t1, t2))
                sum_derivs = 0
                for (j, member2) in enumerate(rec.members)
                    sum_derivs += rec.dg(rec, member.ID, j, t1, t2)
                end
                lyap_condition = (sum_derivs >= 0)
                println("lyap ", lyap_condition)
            end
        end
    end
end

model, objective, rec_centralised, var_members  = optimal_centralised(rec)
println("V(rec)", rec.V(rec))
println("V(rec)", rec.V(rec))
println("Is it NE ", NE_check(rec_centralised))
println(" Summary ", rec_desglose(rec))

iter = length(rec_timeseries)
summary = []
g_members = zeros(iter, N)
P_members = zeros(iter, N)
for j in 1:iter
    data = zeros(T, iter)
    data_gi = zeros(length(rec.members))
    push!(summary, err[min(length(err),j)])
    for (i,member) in enumerate(rec_timeseries[j].members)
        data[:,j] = member.flex_load
        #println("flex load ", i, " iter ", j, " ", member.flex_load)
        data_gi[i] = rec_timeseries[j].g(rec_timeseries[j], member)
    end
    g_members[j,:] = data_gi
    #CSV.write("ResultsIter"*string(j)*".csv", DataFrame(data, :auto), header = false)
end

k = 2
base = 10
for (i,member) in enumerate(rec.members)
    data_Pi = zeros(iter, k)
    for j in 1:iter
        data_Pi[j,:] = rec_timeseries[j].members[i].flex_load[base:base+k-1]
    end
    CSV.write("results/RP"*string(i)*".csv", DataFrame(data_Pi, :auto), header = false)
end

k = 3
base = 10
for (i,member) in enumerate(rec.members)
    data_Pi = zeros(iter, k)
    datatotal_Pi = zeros(iter, T)
    for j in 1:iter
        data_Pi[j,:] = rec_timeseries[j].members[i].flex_load[base:base+k-1]
        datatotal_Pi[j,:] = rec_timeseries[j].members[i].flex_load
    end
    CSV.write("results/3DRP"*string(i)*".csv", DataFrame(data_Pi, :auto), header = false)
    CSV.write("results/RPtotal"*string(i)*".csv", DataFrame(datatotal_Pi, :auto), header = false)
end


data_last = zeros(N,T)

for (i, member) in enumerate(rec.members)
    data_last[i, :] = member.flex_load
end

CSV.write("results/P_last.csv", DataFrame(data_last, :auto), header = false)
CSV.write("results/summary.csv", DataFrame(summary), header = ["error", "revenue", "V(x_k)", "is_NE" , "is_local_maxima"])
CSV.write("results/g_members.csv", DataFrame(g_members, :auto), header = false)



