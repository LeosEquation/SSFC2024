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
function Derivative_IFT(f::Function, x::Float64, p::Float64)
    t = Taylor1(1)
    return -derivative(f(x,p+t))(0.0)/derivative(f(x+t,p))(0.0)
end

#-

"""
    Derivative_IFT(f::Function, x::Float64, p::Vector{Float64},indice::Int64)

Devuelve la derivada `x_p::Float64` de la función implícita de `f`
evaluada en `x` y `p[indice]`.
"""
function Derivative_IFT(f::Function, x::Float64, p::Vector{Float64},indice::Int64)
    t = Taylor1(1)
    T = [i == indice ? t : Taylor1(0) for i in 1:length(p)]
    return -derivative(f(x,p+T))(0.0)/derivative(f(x+t,p))(0.0)
end

#-

"""
    Derivative_IFT(f!::Function, x::Vector{Float64}, p::Float64)

Devuelve la derivada `x_p::Vector{Float64}` de la función implícita del sistema
de ecuaciones diferenciales asicado a `f!` evaluada en `x` y `p`.
"""
function Derivative_IFT(f!::Function, x::Vector{Float64}, p::Float64)
    J = Jacobian(f!,x,p)
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
function Derivative_IFT(f!::Function, x::Vector{Float64}, p::Vector{Float64}, indice::Int64)
    J = Jacobian(f!,x,p,indice)
    Jx = J[:,1:end-1]
    Jp = J[:,end]
    return -inv(Jx)*Jp
end
