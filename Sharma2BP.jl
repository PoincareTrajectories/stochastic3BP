using DifferentialEquations
using Dates

using Plots; #plotly() # Using the Plotly backend

#use GR module
gr();

μ = 10.0
u0 = [1.0,0.0,0.01,1.5]
tspan = (0.0,2.0)

# check if inital conditions exceed escape velocity 
vₑ = √(μ/(u0[1])) #escape velocity
v₀ = √(u0[3]^2 + u0[4]^2 * u0[1]^2) #inital velocity
if vₑ < v₀ 
  throw("Inital velocity is greater than escape velocity")
end

function Sharma2bp(du,u,p,t)
  du[1] = u[3]
  du[2] = u[4] # u2 is angle, [u4]=1/time.
  du[3] = u[1]*u[4]^2 - μ/u[1]^2
  du[4] = -2*u[3]*u[4]/u[1]
end

function σ_2bp(du,u,p,t)
  # du[1] = 0.0
  # du[2] = 0.0
  # du[3] = 0.0121*u[1]
  # du[4] = 0.00022/u[1]

  du[1,1] = 0.0
  du[1,2] = 0.0
  du[2,1] = 0.0
  du[2,2] = 0.0
  du[3,1] = 0.0 # 0.0121*u[1]
  du[3,2] = 0.0
  du[4,1] = 0.0
  du[4,2] = 0.0 # 0.00022/u[1]
end

prob_sde_2bp = SDEProblem(Sharma2bp,σ_2bp, u0, tspan, noise_rate_prototype=zeros(4,2))
sol = solve(prob_sde_2bp, EM(), dt=0.00001)
prob_ode_2bp = ODEProblem(Sharma2bp, u0, tspan)
sol_ode = solve(prob_ode_2bp, RK4(), dt= 0.00001)

# plot(sol,vars=(1,2))
# savefig("test.png")
# # png("test")

r     = sol[1, :]
theta = sol[2, :]
vᵣ    = sol[3, :]
ω     = sol[4, :]
V = sqrt.(sol[3,:].^2 + (sol[4,:].*sol[1,:]).^2)
energy = V.^2 / 2 - μ./r
energy_plot = plot(energy, title="Energy")

x = r.*cos.(theta)
y = r.*sin.(theta)
path_plot = plot(x, y, title ="`SDE Orbit")

r     = sol_ode[1, :]
theta = sol_ode[2, :]
vᵣ    = sol_ode[3, :]
ω     = sol_ode[4, :]
x = r.*cos.(theta)
y = r.*sin.(theta)
path_ode_plot = plot(x, y, title ="ODE Orbit")

#angular momentum
h = r.^2 .*ω

h_plot = plot(h, title="Angular momentum h")

plot(path_plot, path_ode_plot, energy_plot, h_plot, layout=(2,2))
date_string = Dates.format(now(), "YYYY_mm_dd-HH_MM")
u0_string =   replace(string(u0), ['[', ']', ',']=> "")
u0_string = replace(u0_string, [' '] => "_")
savefig("path_energy_h_" * u0_string * date_string * ".png")
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