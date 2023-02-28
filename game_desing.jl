using XLSX
using JuMP
using Gurobi
using CSV
using DataFrames
using Pkg
using Plots, PlotlyJS
using LinearAlgebra
using Roots
using Printf
include("User.jl")
include("bargaining.jl")
include("proximal_gradient_descent.jl")
include("import_data.jl")

Base.show(io::IO, f::Float64) = @printf(io, "%1.2f", f)

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

function best_response_dynamics(rec, type ="proportional", upload = true, learning_rate = 1)
    n_iter = 30
    error = Vector{Float64}(undef, n_iter)
    for i in 1:n_iter
        rec_t = deepcopy(rec)
        new_members = []
        for member in rec.members
            saved_load = member.flex_load
            member.flex_load = 0*member.flex_load
            Set_total_load(rec)
            Set_total_power(rec)
            model, var_member = best_response_i(member, rec, type)
            if upload
                member.flex_load = learning_rate*value.(var_member.P) + (1-learning_rate)*saved_load
            else
                new_member = deepcopy(member)
                new_member.flex_load = learning_rate*value.(var_member.P) + (1-learning_rate)*saved_load
                push!(new_members, new_member)
            end
        end

        if !upload
            rec.members = new_members
        end 
        error[i] = rec_error(rec_t, rec)
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
    lambda_sell = rec.lambda_pun
    lambda_buy = member.lambda_tarif
    println()
    if member.BESS.SOC_max != 0
        println("STORAGE")
    end
    println( "obj value ", objective_value(model) )
    println( "model is ", termination_status(model) )
    println( "Total Load ", rec.load_virtual )
    println( "Total power ", rec.power_virtual )
    P_plus = max.(value.(var_member.P), 0)
    P_minus = max.(-value.(var_member.P), 0)
    println( "costs ", sum(-lambda_buy.*P_plus + lambda_sell.*P_minus) )
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
Pmax = 8
print(lambda_pun)
for i in 1:3
    if i <= 2
        bess = empty_BESS()
        new_member = Member(i, 3*ones(T), zeros(T), 5*Pmax, Pmax, lambda_pun, bess)
    else
        bess = BESS(0, 10, 10)
        new_member = Member(i, zeros(T), zeros(T), Pmax, Pmax, lambda_pun, bess)
    end
    push!(members, new_member)
end

power = 20*rand(Float64, T)
for i in 1:length(power)
    power[i] = min(i, 24-i)
    if (i <= 6) || (i>=18)
        power[i] = 0
    end
end

println("power", power)
lambda_prem =   1
load = 10*rand(Float64, T)
rec = REC(members, power, power, load, load, lambda_pun, lambda_prem)
Set_total_load(rec)
Set_total_power(rec)
println("load", rec.load_virtual)
println("P_res", rec.power_virtual)

#BESTRESPONSEDYNAMICS
type = "proportional"
type = "no"

error = best_response_dynamics(rec, type)
println("error ", error)

display_rec(rec)
display_rec_load(rec)

"""model, rec_centralised = optimal_centralised(rec)
Set_total_load(rec_centralised)
Set_total_power(rec_centralised)
PoA = rec_revenue(rec) - rec_revenue(rec_centralised)
println("PoA ", PoA)
"""
