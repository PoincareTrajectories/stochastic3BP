using Distributions, Random
export solve

function euler_maruyama(μ, σ, t, u0, dt, ΔW, du)
    p = 0.0 #TODO clarify what is p
    return u0 .+ μ(du, u0, p , t) .*dt .+ σ(du, u0, p, t) .* ΔW
end

function solve(μ, σ, u0, tspan; dt=0.01)
    @assert tspan[1]< tspan[2]
    sol = Vector{Any}()
    tvec = Vector{Float64}()
    push!(sol, u0)
    t = tspan[1]
    push!(tvec, t)   
    d =  Normal(0.0, √dt)
    while t < tspan[2]
        u1 = euler_maruyama(μ, σ, t,last(sol),dt, rand(d), 0.0)  
        push!(sol, u1)
        t += dt
        push!(tvec, t)   
    end
    return sol, tvec
end
