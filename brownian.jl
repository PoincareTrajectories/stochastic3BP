# using DifferentialEquations
using Dates

using Plots; 

include("solvers.jl")
#use GR module
gr();

u0 = 0.0
tspan = (0.0,10.0)
function μ(du, u, p , t)
    return sin(t)
end

function σ(du,u,p,t)
    return 1
end

sol, t = solve(μ, σ, u0, tspan, dt = 0.1)
sol2, _ = solve(μ, σ, u0, tspan, dt = 0.1)
sol3, _ = solve(μ, σ, u0, tspan, dt = 0.1)
plot(t, [sol, sol2, sol3])
date_string = Dates.format(now(), "YYYY_mm_dd-HH_MM")
savefig("solutions" * date_string)
