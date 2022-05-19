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

# sol, t = solve(μ, σ, u0, tspan, dt = 0.1)
# sol2, _ = solve(μ, σ, u0, tspan, dt = 0.1)
# sol3, _ = solve(μ, σ, u0, tspan, dt = 0.1)
# sol4, _ = solve(μ, σ, u0, tspan, dt = 0.1)

sols = Vector{Vector{Number}}()
ts = Vector{Any}()
for i in 1:100
    sol, t = solve(μ, σ, u0, tspan, dt = 0.1)
    push!(sols, sol)
    global ts = t
end

plot(ts, sols)
date_string = Dates.format(now(), "YYYY_mm_dd-HH_MM")
savefig("MC-solutions" * date_string)
