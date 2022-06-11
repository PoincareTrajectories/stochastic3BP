using DifferentialEquations
using Dates

using Plots; #plotly() # Using the Plotly backend
using Statistics
using Distributions

#use GR module
gr();


μ = 3.0;
u0 = [1.0,0.0,0.03,1.7];
T = 2.0 ;
tspan = (0.0,T) ;
t_array = range(0,stop=T,length=100) ;

# check if inital conditions exceed escape velocity 
vₑ = √(μ/(u0[1])) #escape velocity
v₀ = √(u0[3]^2 + u0[4]^2 * u0[1]^2) #inital velocity
if vₑ < v₀ 
  throw("Inital velocity is greater than escape velocity")
end

function Sharma2bp(du,u,p,t)
  println(t)
  du[1] = u[3]
  du[2] = u[4] # u2 is angle, [u4]=1/time.
  du[3] = u[1]*u[4]^2 - μ/u[1]^2
  du[4] = -2*u[3]*u[4]/u[1] 
end

prob = ODEProblem(Sharma2bp,u0,tspan)
sol_ode = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
# sol_ode = solve(prob_ode_2bp, RK4(), dt=Δt)
# sol_ode = sol_ode(t_array)

r = sol_ode[1, :]
# print(r)
theta = sol_ode[2, :]
V = sqrt.(sol_ode[3,:].^2 + (sol_ode[4,:].*sol_ode[1,:]).^2)
energy = V.^2 / 2 - μ./r
energy_ode = plot(sol_ode.t, energy, legend=false)
# savefig("deterministic_energy.png")

x = r.*cos.(theta)
y = r.*sin.(theta)
plot(x, y, ylabel="y", xlabel="x", legend=false)
# savefig("deterministic_traj.png")

function σ_2bp(du,u,p,t)
  # du[1] = 0.0
  # du[2] = 0.0
  # du[3] = 0.0  # 0.0121*u[1]
  # du[4] = 0.0  # 0.00022/u[1]

  du[1,1] = 0.0
  du[1,2] = 0.0
  du[2,1] = 0.0
  du[2,2] = 0.0
  du[3,1] = 0.0121*u[1]
  du[3,2] = 0.0
  du[4,1] = 0.0
  du[4,2] = 0.00022/u[1]
end

Δt = 1e-4

prob_sde_2bp = SDEProblem(Sharma2bp,σ_2bp, u0, tspan, noise_rate_prototype=zeros(4,2))
sol = solve(prob_sde_2bp, EM(), dt=Δt)
# sol = sol(t_array)

# plot(sol,vars=(1,2))
# savefig("test.png")
# # png("test")

r     = sol[1, :]
theta = sol[2, :]
vᵣ    = sol[3, :]
ω     = sol[4, :]
V = sqrt.(sol[3,:].^2 + (sol[4,:].*sol[1,:]).^2)
energy = V.^2 / 2 - μ./r
energy_plot = plot(sol.t, energy, ylabel="E(t)", legend=false) #, title="Energy")

x = r.*cos.(theta)
y = r.*sin.(theta)
path_plot = plot(x, y, ylabel="y", xlabel="x", title ="Stochastic", legend=false)

#angular momentum
h = r.^2 .*ω

h_plot = plot(sol.t, h, xlabel="t", ylabel="h(t)", legend=false) #, title="Angular momentum h")

r     = sol_ode[1, :]
theta = sol_ode[2, :]
vᵣ    = sol_ode[3, :]
ω     = sol_ode[4, :]
x = r.*cos.(theta)
y = r.*sin.(theta)
path_ode_plot = plot(x, y, ylabel="y", xlabel="x", title ="Deterministic", legend=false)

#angular momentum
h = r.^2 .*ω

h_ode = plot(sol_ode.t, h, xlabel="t", legend=false)


plot(path_ode_plot, path_plot, energy_ode, energy_plot, h_ode, h_plot, layout=(3,2))
date_string = Dates.format(now(), "YYYY_mm_dd-HH_MM")
u0_string =   replace(string(u0), ['[', ']', ',']=> "")
u0_string = replace(u0_string, [' '] => "_")
savefig("path_energy_h_" * u0_string * date_string * ".png")

nsamples = 2

energy_samples = []
h_samples = []

for i in range(1, stop=nsamples)
  local_sol = solve(prob_sde_2bp, EM(), dt=Δt)
  r     = local_sol[1, :]
  theta = local_sol[2, :]
  vᵣ    = local_sol[3, :]
  ω     = local_sol[4, :]
  V = sqrt.(local_sol[3,:].^2 + (local_sol[4,:].*local_sol[1,:]).^2)
  energy = V.^2 / 2 - μ./r

  h = r.^2 .*ω

  global energy_samples = vcat(energy_samples, energy)
  global h_samples = vcat(h_samples, h)
end

energy_mean = mean(energy_samples, dims=2)
energy_std = std(energy_samples, dims=2)

h_mean = mean(h_samples, dims=2)
h_std = std(h_samples, dims=2)

test = energy_samples
# test = skewness(energy_samples[:, 2])
# print(test)

# n = length(x)
# gif_range = range(1, stop = 150)


# anim = @animate for l in gif_range
#     print(l/150)
#     print(' ')
#     plot(x[1:l], y[1:l])
# end 

# gif(anim, "traj.gif", fps = 50)

# vr = sol[3, :]
# vt = r.*sol[4, :]
# V = sqrt.(vr.^2+vt.^2)
# a = r./(2 .- r.*V.^2 ./ μ)
# plot(a)
# savefig("semimajoraxis.png")

################################################################
# calculating the Hamiltonian directly
################################################################

# sigma_r = 0.0121
# sigma_theta = 0.00022
# function μ_H(du,u,p,t)
#   du[1] = 0.5 *  u[1] ^ 2 * sigma_theta ^2 
# end

# function σ_H(du,u,p,t)
#   du[1,1] = u[3] * sigma_r
#   du[1,2] = u[4] *(u[1] ^ 2)*sigma_theta 
# end

# prob_sde_2bp_hamiltonian = SDEProblem(μ_H, σ_H, u0, tspan, noise_rate_prototype=zeros(1,2))
# sol_hamiltonian = solve(prob_sde_2bp_hamiltonian, EM(), dt=0.001)
# plot(sol_hamiltonian, title="Hamiltonian")

# savefig("hamiltonian" * u0_string * date_string * ".png")

