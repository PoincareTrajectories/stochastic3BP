using DifferentialEquations

μ = 10.0
function Sharma2bp(du,u,p,t)
  du[1] = u[3]
  du[2] = u[4] # u2 is angle, [u4]=1/time.
  du[3] = u[1]*u[4]^2 - μ/u[1]^2
  du[4] = -2*u[3]*u[4]/u[1]
end

function σ_2bp(du,u,p,t)
  du[1] = 0.0
  du[2] = 0.0
  du[3] = 0.0121*u[1]
  du[4] = 0.00022/u[1]
end

prob_sde_2bp = SDEProblem(Sharma2bp,σ_2bp,[1.0,0.0,0.01,1.1],(0.0,1500.0))
sol = solve(prob_sde_2bp)

using Plots; plotly() # Using the Plotly backend

#use GR module
gr();

plot(sol,vars=(1,2))
savefig("test.png")
# png("test")

r = sol[1, :]
theta = sol[2, :]

x = r.*cos.(theta)
y = r.*sin.(theta)
plot(x, y)
savefig("traj.png")

n = length(x)
gif_range = range(1, stop = 150)

anim = @animate for l in gif_range
    print(l/150)
    print(' ')
    plot(x[1:800*l], y[1:800*l])
end 

gif(anim, "traj.gif", fps = 50)

vr = sol[3, :]
vt = r.*sol[4, :]
V = sqrt.(vr.^2+vt.^2)
a = r./(2 .- r.*V.^2 ./ μ)
plot(a)
# savefig('semimajoraxis.png')