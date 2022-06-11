using DifferentialEquations
using Dates

using Plots;

#use GR module
gr();

μ = 3.036e-6;
TU = 3.147e7;
LU = 1;

x_L4 = 0.5 - μ;
y_L4 = sqrt(3)/2;
y_L5 = - y_L4

u0 = [x_L4,y_L4,0.0001,0.01]; #0.004

tspan = (0.0,3000.0)
Δt = 1e-5

function r1(u)
return sqrt((u[1] + μ)^2 + u[2]^2)
end

function r2(u)
    return sqrt((u[1] -1 + μ)^2 + u[2]^2)
end

function S3BP(du,u,p,t)
    dΩx = u[1] - (1 - μ)*(u[1] + μ)/r1(u)^3 - μ*(u[1]-1+μ)/r2(u)^3
    dΩy = u[2] - u[2]*(1-μ)/r1(u)^3 - μ*u[2]/r2(u)^3

    du[1] = u[3]
    du[2] = u[4]
    du[3] = dΩx + 2*u[4]
    du[4] = dΩy - 2*u[3]
  end

prob = ODEProblem(S3BP,u0,tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
# sol = solve(prob, Euler(), dt= Δt)


plot(sol, vars=(1,2), labels="Deterministic Trajectory")
scatter!([(-μ, 0)], markersize=10, markercolor=:yellow, labels="Sun") 
scatter!([(1-μ, 0)], markersize=5, markercolor=:blue, labels="Earth")
scatter!([(x_L4, y_L4)], markersize=3, markercolor=:black, labels="L4")
scatter!([(x_L4, y_L5)], markersize=3, markercolor=:black, labels="L5")
savefig("Deterministic_Trajectory.png")

x = sol[1, :]
y = sol[2, :]
r_1 = sqrt.((x .+ μ).^2 + y.^2)
r_2 = sqrt.((x .-1 .+ μ).^2 + y.^2)
V2 = sol[3, :].^2 .+ sol[4, :].^2

Jacobi = x.^2 .+ y.^2 .+ 2 .* (1 - μ)./r_1 .+ 2 .* μ ./ r_2 .- V2
plot(Jacobi, labels="Deterministic Jacobi")
savefig("Deterministic_Jacobi.png")


σ_r = 1e-4 # 0.0121 * TU^(-3/2);
σ_theta = 1e-4 # 0.00022 * LU * TU^(-3/2);

function σ_3bp(du,u,p,t)
    r = sqrt(u[1]^2 + u[2]^2)
    du[1,1] = 0.0
    du[1,2] = 0.0
    du[2,1] = 0.0
    du[2,2] = 0.0
    du[3,1] =  σ_r*u[1]/r
    du[3,2] = -σ_theta*u[2]
    du[4,1] =  σ_r*u[2]/r
    du[4,2] =  σ_theta*u[1]
end

prob_sde_2bp = SDEProblem(S3BP,σ_3bp, u0, tspan, noise_rate_prototype=zeros(4,2))
sol = solve(prob_sde_2bp, EM(), dt=Δt)
# sol = solve(prob_sde_2bp, SRA3())

plot(sol, vars=(1,2), labels="Stochastic Trajectory")
scatter!([(-μ, 0)], markersize=10, markercolor=:yellow, labels="Sun")
scatter!([(1-μ, 0)], markersize=5, markercolor=:blue, labels="Earth")
scatter!([(x_L4, y_L4)], markersize=3, markercolor=:black, labels="L4")
scatter!([(x_L4, y_L5)], markersize=3, markercolor=:black, labels="L5")
savefig("Stochastic_Trajectory.png")


x = sol[1, :]
y = sol[2, :]
r_1 = sqrt.((x .+ μ).^2 + y.^2)
r_2 = sqrt.((x .-1 .+ μ).^2 + y.^2)
V2 = sol[3, :].^2 .+ sol[4, :].^2

Jacobi = x.^2 .+ y.^2 .+ 2 .* (1 - μ)./r_1 .+ 2 .* μ ./ r_2 .- V2
plot(Jacobi, labels="Stochastic Jacobi")
savefig("Stochastic_Jacobi.png")
