using XLSX
using JuMP
using Gurobi, Ipopt, SCIP
using MathOptInterface
using CSV
using DataFrames
using Pkg
using Plots, PlotlyJS
using LinearAlgebra
using Roots
using Convex
include("User.jl")
include("proximal_gradient_descent.jl")
T= 24

#MINLP branch

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

function optimal_response(rec::REC, member::Member, model=0, type = "v2")

    optimizer = Ipopt.Optimizer
    if isa(model, Int)       
        model = Model(optimizer)
    end
#    set_optimizer_attributes(model, "tol" => 1e-4)
    lambda_pun = rec.lambda_pun
    lambda_prem = rec.lambda_prem

    #WE DEFINE THE VARIABLES
    #The BESS CANT BE CHARGED from the grid
    @variable(model, -member.P_max <= P[1:T] <= member.P_max)
    @variable(model, 0 <= P_plus[1:T] <= member.P_max)
    @variable(model, -member.P_max <= P_minus[1:T] <= 0)

    @variable(model, SE[1:T] )
    SE_without_i = min(rec.power_virtual, rec.load_virtual)
    @variable(model, 0 <= g_plus[1:T] )
    @variable(model, 0 <= g_minus[1:T] )

    @constraint(model, [i=1:T], P[i] == P_plus[i] + P_minus[i])
    @constraint(model, [i=1:T], P[i] <= P_plus[i])
    @constraint(model, [i=1:T], P[i] >= P_minus[i])
    @NLconstraint(model, [i=1:T], P_plus[i]*P_minus[i] >= 0)

    @constraint(model, [i=1:T], SE[i] <= P_plus[i] + rec.load_virtual[i])
    @constraint(model, [i=1:T], SE[i] <= -P_minus[i] + rec.power_virtual[i] )

    #WE ADD CONSTRAINTS RELATIVE TO THE FLEXIBILITY
    P_flexmax = min(member.P_max, member.flex_load_max)
    @variable(model, -P_flexmax <= P_flex[1:T] <= P_flexmax)
    @constraint(model, [i = 1:T], member.P_fix[i] + P_flex[i] >= 0)
    @constraint(model, sum(P_flex) == 0)
    
    #WE ADD CONSTRAINTS RELATIVE TO THE BESS
    if member.BESS.P_max != 0
        #eta = member.BESS.eta 
        eta = 0.9 
        @variable(model, -member.BESS.P_max <= P_bess[1:T] <= member.BESS.P_max)
        @variable(model, 0 <= SOC[1:T] <= member.BESS.SOC_max)
        @variable(model, 0 <= P_bess_plus[1:T] <= member.BESS.P_max)
        @variable(model, -member.BESS.P_max <= P_bess_minus[1:T] <= 0 )
        @NLconstraint(model, [i=1:T], P_bess_plus[i]*P_bess_minus[i] >= 0)


        #@constraint(model, [i=1:T], P_bess_plus[i] <= member.BESS.P_max*u_bess[i])
        #@constraint(model, [i=1:T], (-1+u_bess[i])*member.BESS.P_max <= P_bess_minus[i] )

        @constraint(model, [i=1:T], P_bess[i] == P_bess_plus[i] + P_bess_minus[i])
        @constraint(model, [i=1:T], P_bess[i] <= P_bess_plus[i])
        @constraint(model, [i=1:T], P_bess[i] >= P_bess_minus[i])

        @constraint(model, [i=1:(T-1)], SOC[i+1] == SOC[i] + eta*P_bess_plus[i] + P_bess_minus[i]/eta)
        @constraint(model, SOC[1] == member.BESS.SOC_max)
        @constraint(model, SOC[T] >= 0.8*SOC[1])

    end

    if P_flexmax != 0 && member.BESS.P_max != 0
        @constraint(model, [i=1:T], P[i] == P_bess[i] + P_flex[i] + member.P_fix[i])
    elseif P_flexmax != 0 
        @constraint(model, [i=1:T], P[i] == P_flex[i] + member.P_fix[i])
    else
        @constraint(model, [i=1:T], P[i] == P_bess[i] + member.P_fix[i])
    end

    #WE DEFINE OBJECTIVE CONSTRAINTS
    eps =  0.01


    for i in 1:T    
        @constraint(model, g_plus[i] <= P_plus[i] )
        @constraint(model, g_minus[i] <= -P_minus[i] )

        SE_contribution  = SE
        if type =="v3"
            SE_contribution = SE - SE_without_i
        end
        if rec.load_virtual[i] > eps
            @NLconstraint(model, g_plus[i]*(P_plus[i] + rec.load_virtual[i]) <= ( SE_contribution[i] )*P_plus[i]) 
        else
            @NLconstraint(model, g_plus[i] <= (SE_contribution[i]) )
        end


        SE_contribution  = SE
        if type =="v3"
            SE_contribution = SE - SE_without_i
        end
        if rec.power_virtual[i] > eps
            @NLconstraint(model, g_minus[i]*(-P_minus[i] + rec.power_virtual[i] ) <=  ( SE_contribution[i] )*(-P_minus[i])) 
        else
            @constraint(model, g_minus[i] <= ( SE_contribution[i] ) )
        end
        
    end

    SE_incentive =  (sum(g_plus) + sum(g_minus))*lambda_prem

    cost = 0
    if member.flex_load_max != 0
        cost += sum( lambda_pun.*P_flex ) 
    end
    if member.BESS.P_max !=0
        cost += sum( lambda_pun.*P_bess )
    end 

    regularizer = sum(P[i]^2 for i=1:T)
    regularizer = 0
    @NLobjective(model, Max, SE_incentive - cost + regularizer)
    var_member = Variable_Member(member, P, P_plus, P_minus, g_plus, g_minus)
    if member.BESS.P_max !=0
        var_member = Variable_Member(member, P, P_plus, P_minus, g_plus, g_minus)
    end

    return (model, var_member)
end


function optimal_centralised(rec::REC)
    model = Model(Gurobi.Optimizer)
    var_members = []

    #WE DEFINE THE VARIABLES
    #The BESS CANT BE CHARGED from the grid
    @variable(model, 0 <= P_shared[1:T] )

    for member in rec.members
        P =       @variable(model, [1:T], lower_bound=0, upper_bound=member.BESS.P_max)
        P_flex  = @variable(model, [1:T], lower_bound=-member.BESS.P_max, upper_bound=member.BESS.P_max)

        @constraint(model, [i=1:T], P_flex[i] + member.P_fix[i] >= 0)
        @constraint(model, sum(P_flex) == 0)

        if member.BESS.P_max != 0
            eta = 0.9
            P_bess = @variable(model, [1:T], lower_bound=-member.BESS.P_max, upper_bound=member.BESS.P_max)
            P_bess_plus = @variable(model, [1:T], lower_bound=0, upper_bound=member.BESS.P_max)
            P_bess_minus = @variable(model, [1:T], lower_bound=0, upper_bound=member.BESS.P_max)
            SOC = @variable(model, [1:T], lower_bound=0, upper_bound=member.BESS.SOC_max)


            @constraint(model, [i=1:T], P_bess[i] == P_bess_plus[i] + P_bess_minus[i])
            @constraint(model, [i=1:T], P_bess[i] <= P_bess_plus[i])
            @constraint(model, [i=1:T], -P_bess[i] <= P_bess_minus[i])
            @constraint(model, [i=1:(T-1)], SOC[i+1] == SOC[i] + eta*P_bess_plus[i] - P_bess_minus[i]/eta)
            @constraint(model, SOC[1] == member.BESS.SOC_0)
    
            @constraint(model, [i=1:T], P_bess[i] == P_bess_plus[i] + P_bess_minus[i])
            @constraint(model, [i=1:T], P_bess[i] + P_flex[i] + member.P_fix[i] == P[i])
        else
            @constraint(model, [i=1:T], P_flex[i] + member.P_fix[i] == P[i])
        end

        var_member = Variable_Member(member, P, P, P, 0, 0)
        push!(var_members, var_member)
    end

    expr_plus = 0
    expr_minus = 0
    for var_member in var_members
        expr_plus = expr_plus .+ var_member.P_plus
        expr_minus = expr_plus .+ var_member.P_minus
    end

    @constraint(model, [i=1:T], P_shared[i] <= rec.power_virtual[i] + expr_minus[i] )
    @constraint(model, [i=1:T], P_shared[i] <= expr_plus[i])

    costs = 0
    for var_member in var_members
        costs = costs - sum(rec.lambda_pun.*var_member.P) 
    end
    @objective(model, Max, rec.lambda_prem*sum(P_shared) - costs)

    var_rec = Var_REC(rec, var_members, P_shared , P_shared) 
    return model, var_rec
end


function disagreement_payout(rec::REC, member::Member)
    rec_member = REC([member], 0*range(T), member.flex_load, lambda_pun, 0)
    result_member = optimal_centralised(rec_member)

    new_members = []
    for alt_member in rec.members
        if alt_member != member
            push!(new_members, alt_member)
        end
    end

    rec_minus_member = REC([member], rec.power_virtual, rec.load_virtual-member.flex_load, lambda_pun, 0)
    result_rec = optimal_centralised(rec_minus_member)

    return result_rec, result_member
end

function wardrop_eq(rec::REC)
    A, b =  build_coupled_constraints(rec)
    extended_rec = Extended_REC(A, b, rec)

    #Initialization
    tau = 1
    M = size(A)[1]
    dual_var = Vector{Float64}(undef, M)

    K = 10
    H = 10
    k = 1
    h = 1
    eps1 = 1
    eps2 = 1


    while eps1 > 0.01 && k < K
        eps2=1
        while eps2 > 0.01 && h < H

            #x_i_h+1 = BestResponse(X_i_h, lambda_h)
            #Inner loop until xonvergence
            for member in rec.members
                model, var_member = optimal_response(rec, member)
                original_objective = objective_function(model)

                x = stacked_var(model, var_member)
                N = size(A)[2]/len(rec.members)
                start = N*(member.ID - 1) +1
                A_i = A[:, start:start+N]

                @objective(model, Max, original_objective + dual_var'*A_i*x)
                set_silent(model)
                optimize!(model)
                member.flex_load = value.(var_member.P)
            end

            Set_total_load(rec)
            sigma = rec.load_virtual

            z = (1-1/h)*z + 1/h*(sigma)
            rec.load_virtual = z
            h += 1
        end


        #WE UPDATE THE DUAL VARIABLES
        x = x_tilde
        dual_var = dual_var - tau*(b-A*x)
        dual_var = max.(0, dual_var)
        k = k+1
    end
end

#WE UPDATE THE MODEL IF NECESSARY
function dual_coupled_constraints(model, rec::REC, member::Member, A, b)
    
    M = size(A)[1]
    N = size(A)[2]/len(rec.members)
    start = N*(member.ID - 1) +1
    dual_var = @variable(model, [1:M], lower_bound = 0)
    A_i = A[:, start:start+N]
    obj = objective_function(model)
    @objective(model, Max, obj + dual_var'*A_i)

    return dual_var
end

function build_coupled_constraints(rec::REC)
    #m the number of total variables, k per player
    n = len(rec.members)
    n_dimension = 5
    k = n_dimension*T  
    m = k*len(rec.members)
    n_constraints = 2*m

    v_power = zeros(T, m)
    v_load = zeros(T, m)

    for t in range(T)
        for n_i in range(n)
            v_power[t, n*(t*n_dimension) + n_i*n_dimension + 2] = 1
            v_load[t, n*(t*n_dimension) + n_i*n_dimension + 3] = 1
        end
    end

    #A and b define the constraints Ax >= b
    A = Matrix{Float64}(undef, 0, m)
    b = Vectoe{Float64}(undef, 0)
    for (n_i, member) in enumerate(rec.members)
        for t in range(T)
            A_i = Matrix{Float64}(undef, 0)
            b_i = zeros(4)

            A_sumload = v_load[t]
            A_sumpower = v_power[t]
            A_loadcont = v_load[t]
            A_powercont = v_power[t]

            A_sumload[n*(t*n_dimension) + n_i*n_dimension + 4] = 1
            A_sumpower[n*(t*n_dimension) + n_i*n_dimension + 5] = 1
            A_loadcont[n*(t*n_dimension) + n_i*n_dimension + 6] = 1
            A_powercont[n*(t*n_dimension) + n_i*n_dimension + 7] = 1

            vcat(A, A_sumload)
            vcat(A, A_sumpower)
            vcat(A, A_loadcont)
            vcat(A, A_powercont)

            vcat(A, A_i)
            vcat(b, b_i)
        end
    end

    return A, b
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

function best_response_dynamics(rec, type ="porportional")
    n_iter = 2
    for i in 1:n_iter
        for member in rec.members
            Set_total_load(rec)
            best_response_i(member, rec, type)
        end
    end
    return
end

members = []
T = 24
