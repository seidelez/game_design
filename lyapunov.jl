using XLSX
using JuMP
using Gurobi
using CSV
using DataFrames
using Convex
using Pkg
using Plots
using LinearAlgebra
using IterTools
using Random
include("User.jl")


function update_member(rec::REC, member::Member, learning_rate = 1)
    Set_total_load(rec)
    Set_total_power(rec)
    eta = 0.9
    
    function diff_g(rec::REC, member::Member, x, delta, regularizer = false)
        
        x_plus = max.(0, x)
        x_minus = min.(0, x) 

        Set_total_power(rec)
        Set_total_load(rec)
        P = rec.power_virtual 
        P = rec.power_virtual - max.(0, -member.flex_load) + x_minus
        L = rec.load_virtual 
        L = rec.load_virtual - max.(0, member.flex_load) - x_plus
        SE = min.(L, P)
        D = L - P
        SE_saved = SE
        SE_wo = min.(L + x_plus, P - x_minus)

        

        res_plus = x_plus.*SE/(x_plus+L)
        res_minus = (-x_minus).*SE/(-x_minus+P)
        for i = 1:T
            if x_plus[i] + L[i] == 0
                res_plus[i] = 0
            end
            if -x_minus[i] + P[i] == 0
                res_minus[i] = 0
            end
        end
        
        lambda_prem = rec.lambda_prem
        if rec.configuration =="only charge"
            rec.alpha = 1
        end
        
        if rec.payout == "proportional"
            res = sum(rec.alpha*res_plus + (1-rec.alpha)*res_minus)*lambda_prem
        elseif rec.payout == "marginal"
            res = rec.alpha*sum(SE - SE_wo)*lambda_prem
        elseif rec.payout == "shared"
            res = rec.alpha*sum(SE_saved)*lambda_prem/N
        elseif rec.payout == "penalty"
            penalty_plus = ( D .>= 0).*(x_plus .> 0) 
            penalty_minus = ( D .<= 0).*(x_minus .< 0) 
            res_plus = x_plus
            res_minus = x_minus
            res_plus -= penalty_plus.*(D) 
            res_minus -= penalty_minus.*(-D) 
            res = sum(rec.alpha*res_plus + (1-rec.alpha)*res_minus)*lambda_prem
        else
            res = 0
        end

        res +=  - sum(x.*rec.lambda_pun) 

        if regularizer 
            res += abs(delta)
        end
        return res
    end

    function f_zeros_t1t2(rec::REC, member::Member, t1 , t2)
        x = deepcopy(member.flex_load)
        function flex(delta)
            x[t1] -= delta
            x[t2] += delta

            res = LinearAlgebra.norm1(x - member.P_fix) - member.flex_load_max

            return res
        end

        return flex
    end

    function deriv_delta_t1t2(rec::REC, member::Member, t1 , t2)
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

    random = true
    for t1 = 1:T
        for t2 = 1:T
            if random 
                t1 = rand(1:T)
                t2 = rand(1:T) 
            end
            if t1 != t2

                ordered = t1 < t2
                t1new = min(t1, t2) 
                t2new = max(t1, t2)
                t1 = t1new
                t2 = t2new 
                saved_load = deepcopy(member.flex_load)
                updated_load = deepcopy(member.flex_load)
                set_SOC(member)

                candidates = Dict()
                member.flex_load = 0*member.flex_load
                Set_total_power(rec)
                Set_total_load(rec)
                power_diff = rec.power_virtual - rec.load_virtual
                member.flex_load = saved_load
                P_t1 =  updated_load[t1]
                P_t2 =  updated_load[t2]
                P_max = member.P_max 
                P_min = 0
                if member.BESS.P_max != 0
                    P_min = -member.BESS.P_max
                end

                Delta_max = min(P_t1 - P_min, P_max - P_t2)
                Delta_min = max(P_t1 - P_max, P_min - P_t2)
                #println( "limits 0.1 ", Delta_min, " ", Delta_max)

                oneway = false
                if ordered && oneway == true
                    Delta_min = max(Delta_min, 0)
                end
                if ordered ==0 && oneway == true
                    Delta_max = min(Delta_max, 0)
                end

                #println( "limits 0.2 ", Delta_min, " ", Delta_max)

                if member.BESS.P_max > 0
                    a = minimum(member.BESS.SOC[(t1+1):(t2)])
                    b = (maximum(member.BESS.SOC[(t1+1):(t2)])-member.BESS.SOC_max)
                    Delta_max = min(Delta_max, a)
                    Delta_min = max(Delta_min, b)
                else
                    f_zeros = f_zeros_t1t2(rec, member, t1, t2)
                    try
                        zeros = find_zeros(f_zeros, Delta_min, Delta_max)
                        for zero in zeros
                            if zero < 0
                                Delta_min = max(Delta_min, zero)
                            else
                                Delta_max = min(Delta_max, zero)
                            end
                        end
                    catch
                        _ = 0
                    end
                end

                x_candidates = [0.0, Delta_max, Delta_min, P_t1, -P_t2,
                power_diff[t1],  -power_diff[t1], -power_diff[t2], power_diff[t2]]

                deriv_delta = deriv_delta_t1t2(rec, member, t1, t2)
                #println( "limits 2 ", Delta_min, " ", Delta_max)
                try
                    zeros = find_zeros(deriv_delta, Delta_min, Delta_max)
                    for zero in zeros
                        if true
                            #push!(x_candidates, zero)
                        end
                    end
                catch
                    zeros = [0]
                end

                filter!(x -> x <= Delta_max, x_candidates)
                filter!(x -> x >= Delta_min, x_candidates)
                sort(x_candidates)

                if member.BESS.P_max != 0
                    for cand in x_candidates
                        if cand != 0
                            println(cand, "YES")
                        end
                    end
                end

                for candidate in x_candidates
                    x_delta = deepcopy(member.flex_load)
                    eta = 1
                    if candidate < 0
                        eta = 1/eta
                    end
                    x_delta[t1] = x_delta[t1] - candidate/eta
                    x_delta[t2] = x_delta[t2] + candidate*eta
                    candidates[candidate] = diff_g(rec, member, x_delta, candidate) 
                end

                delta_opt = 0
                try
                    delta_opt = findmax(candidates)[2]
                catch
                    delta_opt = 0
                end
                eta = 1
                updated_load[t1] = updated_load[t1] - delta_opt/eta
                updated_load[t2] = updated_load[t2] + delta_opt*eta

                
                member.flex_load = learning_rate*updated_load + (1-learning_rate)*saved_load

                """model, var_member = optimal_response(rec, member)
                for t = 1:T
                    if t != t1 && t != t2
                        @constraint(model, var_member.P[t] == saved_load[t])
                    end
                end
                set_silent(model)
                optimize!(model)
                member.flex_load = learning_rate*value.(var_member.P) + (1-learning_rate)*saved_load"""
            end
        end
    end

    P_plus = max.(member.flex_load, 0)
    P_minus = min.(member.flex_load, 0)
    var_member = Variable_Member(member, member.flex_load, P_plus, P_minus, P_plus, P_minus)
    return var_member
end

function power_shift(rec, member, t1, t2)
    U_t1 = member.flex_load[t1]
    U_t2 = member.flex_load[t2]
    P1_max = 0
end


function lyapunov(rec::REC, n_iter = 5)
    error = []
    for n_i = 1:n_iter
        rec_2 = deepcopy(rec) 
        for member in rec.members
            var_member = update_member(rec, member)
            #results_summary(0, rec, member, var_member)
        end
        push!(error, rec_error(rec_2, rec))
    end
    return error
end

function projection(rec::REC, member::Member)
    model = Model(Ipopt.optimizer)
    eta = 0.9

    model[:x_projection] = @variable(model, [i=1:T], -member.P_max <= x_projection[i] <= member.P_max)
    @variable(model, [i=1:T], -member.BESS.SOC_max <= SOC[i] <= member.BESS.SOC_max)

    if member.BESS.P_max == 0
        @constraint(model, [i=1:T], x_projection[i]  - member.P_fix[i] >= 0)
    else
        x_proj_plus = Convex.pos(x_projection)
        x_proj_minus = Convex.neg(x_projection)

        @constraint(model, [i=1:(T-1)], SOC[i+1] == SOC[i] + eta*x_proj_plus[i] + x_proj_minus[i]/eta)
        @constraint(model, SOC[T] >= 0.8*SOC[1])
    end

    @NLobjective(model, Min, norm(x_projection, member.flex_load))

    return model
end

