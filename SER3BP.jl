using DifferentialEquations
function f(du,u,p,f,W)
  den1 = (u[1] + mu)^2 + u[2]^2
  den2 = (u[1] -1 + mu)^2 + u[2]^2
  dΩx = u[1] - (1 - mu)*(u[1] + mu)/den1^(3/2) - mu*(u[1]-1+mu)/den2^(3/2)
  dΩy = u[2] - u[2]*(1-mu)/den1^(3/2) - mu*u[2]/den2^(3/2)

  du[1] = u[3]
  du[2] = u[4]
  du[3] = dΩx/(1+e*(1+0.1*W[1])*cos(f)) + 2*u[4]
  du[4] = dΩy/(1+e*(1+0.1*W[1])*cos(f)) - 2*u[3]
end

e = 0.0#0.0489
mu = 0.0009537

u0 = [0.49-mu;sqrt(3)/2; 0.0; 0.0]
tspan = (0.0,100.0)
prob = RODEProblem(f,u0,tspan)
sol = solve(prob,RandomEM(),dt=1/100)

using Plots; plotly() # Using the Plotly backend
plot(sol, vars=(1,2))
