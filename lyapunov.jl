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
include("User.jl")


function update_member(rec::REC, member::Member)
    Set_total_load(rec)
    Set_total_power(rec)
    
    function g(x)
        P = rec.power_virtual
        L = rec.load_virtual
        x_plus = max.(x, 0)
        x_minus = max.(-x, 0)
        K1 = min.(L + x_plus, P)
        K2 = min.(L , P + x_minus)

        res = x_plus.*K1./(x_plus+L)
        res += x_minus.*K2./(x_minus+P)
        return res
    end

    function subdiff_g(x)
        res = zeros(T)
        eps = 1e-3
        for i in 1:T
            e_i = zeros(T)
            e_i[i] =  1
            res[i] = (sum(g(x+eps*e_i))-sum(g(x-eps*e_i)))/(2*eps)
        end
        return subdiff_g
    end

    model = projection(rec, member)
    optimize!(model)
    member.flex_load = value.(model[:x_projection])

    return model, g(member.flex_load)
end

function projection(rec::REC, member::Member)
    model = Model(Ipopt.optimizer)
    eta = 0.9

    model[:x_projection] = @variable(model, [i=1:T], -member.P_max <= x_projection[i] <= member.P_max)
    @variable(model, [i=1:T], -member.BESS.SOC_max <= SOC[i] <= member.BESS.SOC_max)

    if member.BESS.P_max == 0
        @constraint(model, [i=1:T], x_projection[i]  - member.P_fix[i] >= 0)
    else
        @expression(model, x_proj_plus = Convex.pos(x_projection))
        @expression(model, x_proj_minus = Convex.neg(x_projection))

        @constraint(model, [i=1:(T-1)], SOC[i+1] == SOC[i] + eta*x_proj_plus[i] + x_proj_minus[i]/eta)
        @constraint(model, SOC[T] >= 0.8*SOC[1])
    end

    @NLobjective(model, norm(x_projection, member.flex_load))

    return model
end

