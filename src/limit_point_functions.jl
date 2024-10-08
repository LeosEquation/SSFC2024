# # Funciones de puntos límite
# Este archivo contiene funciones correspondientes al sistema
# de punto límite así como su jacobiano y los valores iniciales

#-

function LP_Function!(F::Vector{Float64}, f::Function, dfx::Float64, x::Float64, p::Union{Float64, Vector{Float64}}, v::Float64, w::Float64, t::Float64)

    F[1] = f(x, p, t)
    F[2] = dfx*v
    F[3] = w*v - 1.0

end

#-

function LP_Jacobian!(J::Matrix{Float64}, dfx::Float64, dfp::Float64, dfxx::Float64, dfxp::Float64, v::Float64, w::Float64)

    J[1,1] = dfx
    J[1,2] = dfp

    J[2,1] = dfxx*v
    J[2,2] = dfxp*v
    J[2,3] = dfx

    J[3,3] = w

end

#-

function LP_Function!(F::Vector{Float64}, dfx::Matrix{Float64}, dx::Vector{Float64}, v::Vector{Float64}, w::Vector{Float64}, n::Int64)

    F[1:n]         .= dx
    F[(n+1):(2*n)] .= dfx*v
    F[2*n + 1]      = w ⋅ v - 1.0

end

#-

function LP_Jacobian!(J::Matrix{Float64}, dfx::Matrix{Float64}, dfp::Vector{Float64}, dfxx::Array{Float64, 3}, dfxp::Matrix{Float64}, v::Vector{Float64}, w::Vector{Float64}, n::Int64)

    ###

    J[1:n, 1:n] .= dfx

    J[1:n, n+1] .= dfp

    #J[1:n, (n+2):(2*n+1)] .= zeros(n,n)

    ####

    J[(n+1):(2*n),1:n]           .= [sum([dfxx[i, k, j]*v[k] for k in 1:n]) for i in 1:n, j in 1:n]

    J[(n+1):(2*n),n+1]           .= [sum([dfxp[i, k]*v[k] for k in 1:n])    for i in 1:n]

    J[(n+1):(2*n),(n+2):(2*n+1)] .= dfx

    ###

    #J[2*n+1, 1:n] .= zeros(n)

    #J[2*n+1, n+1] .= 0.0

    J[2*n+1, (n+2):(2*n+1)] .= w

end


