# # Función Implícita
# Este archivo contiene las funciones necesarias para calcular
# la derivada de la función implícita haciendo uso del Teorema
# de la función Implícita.

#-

"""
    Derivative_IFT(f::Function, x::Float64, p::Float64)

Devuelve la derivada `x_p::Float64` de la función implícita de `f`
evaluada en `x` y `p`.
"""
function Derivative_IFT(f::Function, x::Float64, p::Float64,t::Float64)
    return -derivative(f(x+r,p+s,t+r))(0.0)/derivative(f(x+s,p+r,t+r))(0.0)
end

#-

"""
    Derivative_IFT(f::Function, x::Float64, p::Vector{Float64},indice::Int64)

Devuelve la derivada `x_p::Float64` de la función implícita de `f`
evaluada en `x` y `p[indice]`.
"""
function Derivative_IFT(f::Function, x::Float64, p::Vector{Float64},t::Float64,indice::Int64)
    return -derivative(f(x+r,p+S,t+r))(0.0)/derivative(f(x+s,p .+ r,t+r))(0.0)
end

#-

"""
    Derivative_IFT(f!::Function, x::Vector{Float64}, p::Float64)

Devuelve la derivada `x_p::Vector{Float64}` de la función implícita del sistema
de ecuaciones diferenciales asicado a `f!` evaluada en `x` y `p`.
"""
function Derivative_IFT(f!::Function, x::Vector{Float64}, p::Float64,t::Float64)
    J = Jacobian(f!,x,p,t)
    Jx = J[:,1:end-1]
    Jp = J[:,end]
    return -inv(Jx)*Jp
end

#-

"""
    Derivative_IFT(f!::Function, x::Vector{Float64}, p::Vector{Float64},indice::Int64)

Devuelve la derivada `x_p::Vector{Float64}` de la función implícita del sistema
de ecuaciones diferenciales asicado a `f!` evaluada en `x` y `p[indice]`.
"""
function Derivative_IFT(f!::Function, x::Vector{Float64}, p::Vector{Float64},t::Float64, indice::Int64)
    J = Jacobian(f!,x,p,t,indice)
    Jx = J[:,1:end-1]
    Jp = J[:,end]
    return -inv(Jx)*Jp
end
