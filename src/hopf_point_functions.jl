# # Funciones de puntos de Hopf
# Este archivo contiene funciones correspondientes al sistema
# de Hopf así como su jacobiano y los valores iniciales

#-

#using TaylorSeries, LinearAlgebra

#-

function Hopf_Function!(F::Vector{Float64}, dx::Vector{Float64}, dfx::Matrix{Float64}, v::Vector{Float64}, w::Vector{Float64}, k::Float64, n::Int64)

    F[1:n] .= dx
    F[n + 1:2*n] .= dfx^2 * v + k * v
    F[2*n + 1] = v ⋅ v - 1.0
    F[end] = w ⋅ v 

end

function Hopf_Jacobian!(J::Matrix{Float64}, dfx::Matrix{Float64}, dfp::Vector{Float64}, dfxx::Array{Float64, 3}, dfxp::Matrix{Float64}, v::Vector{Float64}, w::Vector{Float64}, k::Float64, n::Int64)

    J[1:n, 1:n] .= dfx
    J[1:n, n+1] .= dfp
    
    J[n+1:2*n, 1:n] .= [sum([(dfxx[i,l,j]*dfx[l,k] + dfx[i,l]*dfxx[l,k,j])*v[k] for k in 1:n, l in 1:n]) for i in 1:n, j in 1:n]
    J[n+1:2*n, n+1] .= [sum([(dfxp[i,l]*dfx[l,k] + dfx[i,l]*dfxp[l,k])*v[k] for k in 1:n, l in 1:n]) for i in 1:n]
    J[n+1:2*n, n+2:2*n+1] .= dfx^2 + k*I(n)
    J[n+1:2*n, end] .= v

    J[end-1, n+2:2*n+1] .= 2*v

    J[end, n+2:2*n+1]   .= w

end

#-
