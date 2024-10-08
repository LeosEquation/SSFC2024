




function Equilibrium_Function!(F::Vector{Float64}, f::Function, x::Float64, p::Float64, t::Float64, x0::Float64, p0::Float64, xs::Float64, ps::Float64, Δs::Float64)
    
    F[1] = f(x,p,t)
    F[2] = (x - x0) * xs + (p - p0) * ps - Δs

end

function Equilibrium_Function!(F::Vector{Float64}, f::Function, x::Float64, p::Vector{Float64}, t::Float64, x0::Float64, p0::Vector{Float64}, xs::Float64, ps::Float64, Δs::Float64, indice::Int64)
    
    F[1] = f(x,p,t)
    F[2] = (x - x0) * xs + (p[indice] - p0[indice]) * ps - Δs

end

function Equilibrium_Function!(F::Vector{Float64}, f!::Function, x::Vector{Float64}, p::Float64, t::Float64, x0::Vector{Float64}, p0::Float64, xs::Vector{Float64}, ps::Float64, Δs::Float64, dx::Vector{Float64}, n::Int64)
    
    f!(dx, x, p, t)

    F[1:n] .= dx
    F[n+1] = (x - x0) ⋅ xs + (p - p0) * ps - Δs

end

function Equilibrium_Function!(F::Vector{Float64}, f!::Function, x::Vector{Float64}, p::Vector{Float64}, t::Float64, x0::Vector{Float64}, p0::Vector{Float64}, xs::Vector{Float64}, ps::Float64, Δs::Float64, dx::Vector{Float64}, n::Int64, indice::Int64)
    
    f!(dx, x, p, t)

    F[1:n] .= dx
    F[n+1] = (x - x0) ⋅ xs + (p[indice] - p0[indice]) * ps - Δs

end



function Equilibrium_Jacobian!(J::Matrix{Float64}, f::Function, x::Float64, p::Float64, t::Float64, xs::Float64, ps::Float64, s::Taylor1{Float64}, r::Taylor1{Float64})

    J[1,1] = differentiate(f(x+s,p+r,t+r))(0.0)
    J[2,1] = xs
    J[1,2] = differentiate(f(x+r,p+s,t+r))(0.0)
    J[2,2] = ps

end

function Equilibrium_Jacobian!(J::Matrix{Float64}, f::Function, x::Float64, p::Vector{Float64}, t::Float64, xs::Float64, ps::Float64, s::Taylor1{Float64}, r::Taylor1{Float64}, S::Vector{Taylor1{Float64}})

    J[1,1] = differentiate(f(x + s, p .+ r, t + r))(0.0)
    J[2,1] = xs
    J[1,2] = differentiate(f(x + r, p .+ S, t + r))(0.0)
    J[2,2] = ps

end

function Equilibrium_Jacobian!(J::Matrix{Float64}, f!::Function, x::Vector{Float64}, p::Float64, t::Float64, xs::Vector{Float64}, ps::Float64, dx::Vector{Taylor1{Float64}}, s::Taylor1{Float64}, r::Taylor1{Float64}, S::Matrix{Taylor1{Float64}}, n::Int64)

    for i in 1:n
        f!(dx, x .+ S[:,i], p + r, t + r)
        J[1:n,i] .= differentiate.(dx)(0.0)
    end

    J[end, 1:n] .= xs

    f!(dx, x .+ r, p + s, t + r)

    J[1:n, n+1] .= differentiate.(dx)(0.0)

    J[end, end] = ps

end

function Equilibrium_Jacobian!(J::Matrix{Float64}, f!::Function, x::Vector{Float64}, p::Vector{Float64}, t::Float64, xs::Vector{Float64}, ps::Float64, dx::Vector{Taylor1{Float64}}, s::Taylor1{Float64}, r::Taylor1{Float64}, Sx::Matrix{Taylor1{Float64}}, Sp::Vector{Taylor1{Float64}}, n::Int64)

    for i in 1:n
        f!(dx, x .+ Sx[:,i], p .+ r, t + r)
        J[1:n,i] .= differentiate.(dx)(0.0)
    end

    J[end, 1:n] .= xs

    f!(dx, x .+ r, p .+ Sp, t + r)

    J[1:n, n+1] .= differentiate.(dx)(0.0)

    J[end, end] = ps

end