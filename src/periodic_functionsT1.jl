

#-

function Function_Periodic_Orbits!(F::Vector{Float64}, x::Vector{Float64}, p::Float64, T::Float64, x0::Vector{Float64}, p0::Float64, T0::Float64, xs::Vector{Float64}, ps::Float64, Ts::Float64, Δs::Float64, sol::Vector{Float64}, dx0::Vector{Float64}, n::Int64)
    
    ### \int_{0}^{T} F(x,p) dt

    for i in 1:n

        F[i] = sol[i] - x[i]

    end

    ### Poincare Phase Condition (x(0) - x_{0}(0)) \cdot x'_{0}(0)
    
    F[n+1] = (x - x0) ⋅ dx0

    ### PseudoArc-Length Continuation

    F[n+2] = (x - x0) ⋅ xs + (p - p0)*ps + (T - T0)*Ts - Δs
    
end

#-

function Jacobian_Periodic_Orbits!(J::Matrix{Float64}, f!::Function, p, T::Float64, xs::Vector{Float64}, ps::Float64, Ts::Float64, Φ::Matrix{Float64}, sol::Vector{Float64}, dx::Vector{Float64}, dx0::Vector{Float64}, n::Int64)

    f!(dx, sol, p, T)

    ### Jacobian of \int_{0}^{T} F(x,p) dt

    for i in 1:n

        for j in 1:n

            J[i, j] = Φ[i,j] - ==(i,j)

        end

        J[i, n+1]   = Φ[i, n+1]

        J[i, n+2]   = dx[i]

    end

    ### Jacobian of Poincare Phase Condition (x(0) - x_{0}(0)) \cdot x'_{0}(0)

    for j in 1:n

        J[n+1, j] = dx0[j]

    end

    #   J[n+1, n+1]  = 0.0

    #   j[n+1, n+2]  = 0.0

    ### Jacobian of PseudoArc-Length Continuation

    for j in 1:n

        J[n+2, j]    = xs[j]

    end

        J[n+2, n+1]  = ps

        J[n+2, n+2]  = Ts
    
end

#-