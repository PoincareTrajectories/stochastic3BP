using DifferentialEquations

function Sharma2bp(du,u,p,t)
  du[1] = u[3]
  du[2] = u[4]
  du[3] = u[1]*u[4]^2 - 10.0/u[1]^2
  du[4] = -2*u[3]*u[4]/u[1]
end

function σ_2bp(du,u,p,t)
  du[1] = 0.0
  du[2] = 0.0
  du[3] = 0.0121*u[1]
  du[4] = 0.00022/u[1]
end

prob_sde_2bp = SDEProblem(Sharma2bp,σ_2bp,[1.0,0.0,0.01,1.1],(0.0,150.0))
sol = solve(prob_sde_2bp)

using Plots; plotly() # Using the Plotly backend
plot(sol,vars=(1,2))
