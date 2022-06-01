using Distributions, Random
export solve

function euler_maruyama(μ, σ, t, u0, dt, ΔW, du)
    # @show μ(du, u0, t)
    return u0 + μ(du, u0, t) *dt + σ(du, u0,  t) .* ΔW
end

function ssrkw1(μ, σ, t, u0, dt, ΔW, du)
    println("solving with new approach")
    return du
end


function get_solver(s)
    s == :EM && return euler_maruyama
    s == :ssrkw1 && return ssrkw1
    throw("No such solver implemented")
end


function solve(μ, σ, u0, tspan, s=:EM; dt=0.01, noise_rate_prototype=zeros(1,1))
    @assert tspan[1]< tspan[2]
    solver = get_solver(s)
    r = range(;start=tspan[1],stop=tspan[2], step=dt)
    n = length(r)
    u = zeros(n, length(u0))
    tvec = collect(r)
    u[1,:] = u0
    d =  Normal(0.0, √dt)
    for i in 1:n-1
        u[i+1,:] = solver(μ, σ, tvec[i],u[i,:],dt, rand(d, size(noise_rate_prototype)), Vector{Float64}(undef, length(u0)))  
    end
    return u, tvec
end
