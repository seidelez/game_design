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
    P_max 
end

mutable struct Variable_BESS
    bess
    SoC
    P_max 
end

mutable struct REC
    members
    power_fix
    power_virtual
    load_fix
    load_virtual
    lambda_pun
    lambda_prem
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

function empty_BESS()
   return BESS(0, 0, 0)
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
        power_res .+= max.(-member.flex_load, 0)
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
    for member in rec.members
        load_saved = model.flex_load
        model, var_member = optimal_response(rec, member)
        set_silent(model)
        optimize!(model)
        if norm(load_saved - value.(var_member.P)) > 0.001
            return false
        end
    end
    return true
end

function rec_error(rec1::REC, rec2::REC)
    res = 0
    function g(rec)
        P = rec.power_virtual
        L = rec.load_virtual
        K1 = min.(L , P)
        K2 = min.(L , P)
        SE = min.(L , P)

        res = sum(SE)*rec2.lambda_prem - sum(L.*rec2.lambda_pun)
        return res
    end


    for (i,member) in enumerate(rec1.members)
        v = rec1.members[i].flex_load - rec2.members[i].flex_load
        res +=  LinearAlgebra.norm(v)
    end

    res = g(rec2)
    return res
end

function rec_revene(rec)
    SE = min.(rec.load_virtual, rec.power_virtual)
    revenue = rec.lambda_prem*sum(SE) - sum(rec.load_virtual.*rec.lambda_pun) 
end

