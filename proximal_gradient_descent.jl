using JuMP
#using Gurobi
using CSV
using DataFrames
using Pkg
using Plots, PlotlyJS
using LinearAlgebra
using Roots
using IterTools

struct polytope
    A
    b
end

struct prox_grad_problem
    f
    g
    gamma
end

function prox_operator(polytope)

    function projection(x)
        if polytope.A*x <= polytope.b
            return x
        end

        proj = 0
        return proj
    end

    return projection
end

##inicio con add_elem(n)
##if elem(n) > n add_elem(n-1)
function add_elem(elem, i, n)
    siz = length(elem)
    j = siz-i
    mark = 0
    elem[i] = elem[i] + 1

    if i == 1 && elem[i]>n-j
        add_elem(elem, i+1, n)
        return false
    end

    if elem[i] > n-j
        elem[i+1] = elem[i] + 1

        add_elem(elem, i-1, n)
        return
    end

    return
end 
