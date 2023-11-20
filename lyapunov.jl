using XLSX
using JuMP
#using Gurobi
using CSV
using DataFrames
using Convex
using Pkg
using Plots
using LinearAlgebra
using IterTools
using Random
include("User.jl")


function update_member(rec::REC, member::Member, t1, t2, learning_rate = 1)
    Set_total_load(rec)
    Set_total_power(rec)
    eta = 0.9
    
    function diff_g(rec::REC, member::Member, x)
        updated_member = deepcopy(member)
        updated_member.flex_load = updated_member.flex_load + x

        res = rec.g(rec, member)
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

    saved_load = deepcopy(member.flex_load)
    updated_load = deepcopy(member.flex_load)
    set_SOC(member)
    candidates = Dict()
    member.flex_load = 0*member.flex_load
    Set_total_power(rec)
    Set_total_load(rec)
    power_diff = rec.power_virtual - rec.load_virtual
    member.flex_load = saved_load
        
    Delta = rec.C(rec, member, t1, t2)
    x_candidates = [Delta.a, 0, Delta.b]
    try
        zeros = find_zeros(deriv_delta, Delta_min, Delta_max)
        for zero in zeros
            push!(x_candidates, zero)
        end
    catch
        zeros = [0]
    end
    sort(x_candidates)
    for candidate in x_candidates
        member_candidate = deepcopy(member)
        eta = 1
        if candidate < 0
            eta = 1/eta
        end
        member_candidate.flex_load = member_candidate.flex_load + candidate*v_power(t1,t2)
        candidates[candidate] = rec.g(rec, member_candidate)
    end
    delta_opt = 0
    try
        delta_opt = findmax(candidates)[2]
    catch
        delta_opt = 0
    end
        eta = 1
    
    updated_load = updated_load + v_power(t1, t2)*delta_opt/eta
    
        
    member.flex_load = learning_rate*updated_load + (1-learning_rate)*saved_load
    
    #println("P diff: ", sum(member.flex_load))
    #println("delta opt: ", delta_opt)
    #println("delta opt: ", Delta)
    var_member = 0
    return delta_opt, var_member
end

function power_shift(rec, member, t1, t2)
    U_t1 = member.flex_load[t1]
    U_t2 = member.flex_load[t2]
    P1_max = 0
end

function lyapunov(rec::REC;  n_iter = 10)
    error = []
    if rec.single_update 
        n_iter = 50
    end
    rec_timeseries = [rec]
    single_update = rec.single_update
    for n_i = 1:n_iter
        rec_2 = deepcopy(rec) 
        println("NEW ITER")
        
        if !single_update
            for member in rec.members
                for t1 in 1:T
                    for t2 in 1:T
                        delta_star, var_member = update_member(rec, member, t1, t2)
                    end
                    new_error = rec_error(rec_2, rec)
                    push!(error, new_error)
                    push!(rec_timeseries, rec)
                end
            end
            if equal_rec(rec, rec_2)
                return error, rec_timeseries
            end
        else
            candidate_member = Dict()
            for (i, member) in enumerate(rec.members)
                for t1 in 1:T
                    for t2 in 1:T
                        rec_copy = deepcopy(rec)
                        member_copy = deepcopy(member)
                        delta_star, var_member = update_member(rec_copy, member_copy, t1, t2)
                        candidate_member[(i, t1, t2)] = rec_copy.g(rec_copy, member_copy) - rec.g(rec, member)
                    end
                end
            end

            rec_copy = deepcopy(rec)
            rec_final = deepcopy(rec)
            cardinal_candidates = 1
            max_comb = findmax(candidate_member) 
            max_val = max_comb[1]
            g_ini = rec.g(rec, rec.members[max_comb[2][1]]) 
            candidate_member[max_comb[1]] = 0

            while max_val == max_comb[1] && max_val != 0
                i, t1, t2 = max_comb[2]
                delta_star, var_member = update_member(rec_copy, rec_copy.members[i], t1, t2)
                for (j, member) in enumerate(rec_final.members)
                    member.flex_load += rec_copy.members[j].flex_load
                end
                cardinal_candidates += 1
                candidate_member[max_comb[2]] = 0
                max_comb = findmax(candidate_member) 
                rec_copy = deepcopy(rec)
            end

            for member in rec_final.members
                member.flex_load /= cardinal_candidates
            end

            rec = rec_final
            new_error = rec_error(rec_2, rec)
            push!(error, new_error)
            push!(rec_timeseries, rec)
        end
        equal = equal_rec(rec, rec_2)
    
        if equal !! new_error[1]==0
            println(" equal ")
            break
        else
            continue
        end
    end
    return error, rec_timeseries
end

function fixed_C(rec::REC)
    error = []
    rec_timeseries = [rec]
    for n_i = 1:n_iter
        rec_2 = deepcopy(rec) 
        println("NEW ITER")
        
        if !single_update
            rec_saved = deepcopy(rec)
            for (i, member) in enumerate(rec.members)
                for t1 in 1:T
                    for t2 in 1:T
                        if member.C(rec, member, t1, t2) == rec_saved.C(rec_saved, rec_save.members[i], t1, t2)
                            delta_star, var_member = update_member(rec, member, t1, t2)
                        end
                    end
                end
            end
            push!(error, rec_error(rec_2, rec))
            push!(rec_timeseries, rec)
        end
        equal = equal_rec(rec, rec_2)
    
        if equal
            println(" equal ")
            break
        else
            continue
        end
    end
    return error, rec_timeseries
end


function markov_lyapunov(rec::REC, n_iter = 10)
    error = []
    for n_i = 1:n_iter
        rec_2 = deepcopy(rec) 
        for member in rec.members
            for t1 in 1:T
                for t2 in 1:T
                    var_member = update_member(rec, member, t1, t2)
                end
            end
            #results_summary(0, rec, member, var_member)
        end
        push!(error, rec_error(rec_2, rec))
        equal = equal_rec(rec, rec_2)
    
        if equal
            break
        else
            continue
        end
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

function displaced_lyapunov(rec::REC, n_iter =5)
    error = []
    rec_2 = deepcopy(rec)
    n = 1 
        
    while !equal_rec(rec, rec_2) || n<n_iter
        n += 1
        rec_2 = deepcopy(rec)
        lyapunov(rec)
        candidates = border_direction(rec)
        candidates2 = border_direction_alt1(rec)
        candidates = vcat(candidates, candidates2)

        for direction in candidates
            rec_3 = deepcopy(rec)
            for (i, member) in enumerate(rec_3.members)
                member.flex_load += direction[i, :]
            end

            feasible = true
            feasible = is_REC_feasible(rec_3)
            if feasible && (V(rec_3) > V(rec))
                rec = deepcopy(rec_3)
                println("newrec!!")
                break
            end
        end

        push!(error, rec_error(rec_2, rec))
    end

    return error
end

