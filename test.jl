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
lambda_pun = 0*lambda_pun

lambda_pun[16:24] = lambda_pun[16:24]
lambda_prem = 0.12
lambda_prem = 1*lambda_prem
payout = "marginal"
configuration = "standalone"
alpha = 0.3
load = ones(T)
single_update = false
delta = delta_limited
rec = REC(members2, payout, configuration, power, power, load, load, lambda_pun, lambda_prem, alpha, g_marginal, dg_marginal, delta, V_hybrid, single_update)
Set_total_load(rec)
Set_total_power(rec)
rec_centralised = deepcopy(rec)
rec_initial = deepcopy(rec)
println("load", rec.load_virtual)
println("P_res", rec.power_virtual)


# 
rec_centralised = deepcopy(rec)
model, objective, rec_centralised, var_members = optimal_centralised(rec_centralised)
for (i, member) in enumerate(rec_centralised.members)
    member.flex_load = value.(var_members[i].P)
end
summary_rec(rec_centralised)
println(" objective central", value.(objective))

P = rec.power_virtual
L = rec.load_virtual
println("SE ", min.(L, P))
println("P ", P)
println("L ", L)
println("summary ", rec_error(rec_centralised, rec_centralised))

println("NE check ", NE_check(rec_centralised))

rec2 = rec_centralised
rec2.members = []
for (i,member) in enumerate(rec2.members)
    println("Player: ",i)
    println("initial g", rec2.g(rec2, member))
    println("initial P", member.flex_load)
    println("initial V", rec2.V(rec2))
    model, objective, var_member =  optimal_response(rec2, member)
    saved_P = member.flex_load
    set_silent(model)
    optimize!(model)
    member.flex_load = value.(var_member.P)

end
