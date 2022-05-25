# using DifferentialEquations
using Dates

using Plots; 

include("solvers.jl")
#use GR module
gr();

u0 = [0.0]
tspan = (0.0,10.0)
function μ(du, u, t)
    du[1] = sin(t)
    return du
end

function σ(du,u,t)
    du[1] = 1
    return du
end

# sol, t = solve(μ, σ, u0, tspan, dt = 0.1)
# plot(t, sol)
sols = Vector{Any}()
ts = Vector{Any}()
for i in 1:100
    sol, t = solve(μ, σ, u0, tspan, dt = 0.1)
    push!(sols, sol)
    global ts = t
end

plot(ts, sols)
date_string = Dates.format(now(), "YYYY_mm_dd-HH_MM")
savefig("MC-solutions" * date_string)

################################################################
# running the Lotka-Volterra SDE
################################################################

u0 = [2.5, 45.0]
tspan = (0.0, 10.0)

function mu(du, u, t)
    du[1] = u[1] * ( 1- 0.02 * u[2] )
    du[2] = u[2] * ( 0.01 * u[1] - 1)
    return du
end

function sigma(du, u, t)
    du[1] = 5
    du[2] = 5
    return du
end

lotka_volterra, t = solve(mu, sigma, u0, tspan, :ssrkw1, dt=0.01)
@show size(lotka_volterra)
plot(t, lotka_volterra)
savefig("Lotka Volterra" * date_string)

