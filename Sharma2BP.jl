using DifferentialEquations
using Dates

using Plots; #plotly() # Using the Plotly backend

#use GR module
gr();

μ = 10.0
u0 = [1.0,0.0,0.01,1.1]
tspan = (0.0,1.0)

# check if inital conditions exceed escape velocity 
vₑ = √(μ/(u0[1]))
v₀ = √(u0[3]^2 + u0[4]^2 * u0[1]^2)
if vₑ < v₀ 
  throw("Inital velocity is greater than escape velocity")
end

function Sharma2bp(du,u,p,t)
  du[1] = u[3]
  du[2] = u[4] # u2 is angle, [u4]=1/time.
  du[3] = u[1]*u[4]^2 - μ/u[1]^2
  du[4] = -2*u[3]*u[4]/u[1]
end

# prob = ODEProblem(Sharma2bp,u0,tspan)
# sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# r = sol[1, :]
# theta = sol[2, :]

# x = r.*cos.(theta)
# y = r.*sin.(theta)

# plot(x, y)

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
sol = solve(prob_sde_2bp, EM(), dt=0.000001)

# plot(sol,vars=(1,2))
# savefig("test.png")
# # png("test")

r = sol[1, :]
theta = sol[2, :]
V = sqrt.(sol[3,:].^2 + (sol[4,:].*sol[1,:]).^2)
energy = V.^2 / 2 - μ./r
energy_plot = plot(energy)

x = r.*cos.(theta)
y = r.*sin.(theta)
path_plot = plot(x, y)

plot(path_plot, energy_plot, layoug=(2,1))
date_string = Dates.format(now(), "YYYY_mm_dd-HH_MM")
savefig("path_energy_" * date_string * ".png")
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