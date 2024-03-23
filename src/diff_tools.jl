# # Herramientas para derivación
# Este archivo contiene las funciones necesarias para calcular
# gradientes y jacobianos de las funciónes que representan ecuaciones
# diferenciales.

#-

"""
    Gradient(f::Function,x::Float64,p::Float64)

Devuelve el gradiente `∇f::Vector{Float64}` de la función `f` evaluado
en `x` y `p`.
"""
function Gradient(f::Function,x::Float64,p::Float64)
    t = Taylor1(1)
    return [derivative(f(x+t,p))(0.0),derivative(f(x,p+t))(0.0)]
end

#-

"""
    Gradient(f::Function,x::Float64,p::Vector{Float64},indice::Int64)

Devuelve el gradiente `∇f::Vector{Float64}` de la función `f` evaluado
en `x` y `p[indice]`.
"""
function Gradient(f::Function,x::Float64,p::Vector{Float64},indice::Int64)
    t = Taylor1(1)
    T = [i == indice ? t : Taylor1(0) for i in 1:length(p)]
    return [derivative(f(x+t,p))(0.0); [i == indice ? derivative(f(x,p+T))(0.0) : 0.0 for i in 1:length(p)]]
end

#-

"""
    Gradient(f::Function,x::Vector{Float64},p::Float64)

Devuelve el gradiente `∇||f||::Vector{Float64}` de la norma del sistema
de ecuaciones diferenciales asociado a función `f!` evaluado en `x` y `p`.
"""
function Gradient(f!::Function,x::Vector{Float64},p::Float64)
    t = Taylor1(1)
    s = Taylor1([0.0,0.0],1)
    dx = [s for i in 1:length(x)]
    G = zeros(length(x) + 1)
    for i in 1:length(x)
        f!(dx,x + [i == j ? t : s for j in 1:length(x)],p)
        G[i] = derivative(dx ⋅ dx)(0.0)
    end
    f!(dx,x,p + t)
    G[end] = (dx ⋅ dx)(0.0)
    return G
end

#-

"""
    Gradient(f!::Function,x::Vector{Float64},p::Vector{Float64},indice::Int64)

Devuelve el gradiente `∇||f||::Vector{Float64}` de la norma del sistema
de ecuaciones diferenciales asociado a función `f!` evaluado en `x` y `p[indice]`.
"""
function Gradient(f!::Function,x::Vector{Float64},p::Vector{Float64},indice::Int64)

    t = Taylor1(1)
    s = Taylor1([0.0,0.0],1)
    dx = [s for i in 1:length(x)]
    Gx = zeros(length(x))
    Gp = zeros(length(p))
    T = [i == indice ? t : s for i in 1:length(p)]

    for i in 1:length(x)
        f!(dx,x + [i == j ? t : s for j in 1:length(x)],p)
        Gx[i] = derivative(dx ⋅ dx)(0.0)
    end

    for i in 1:length(p)
        f!(dx,x,p + T)
        if i == indice
            Gp[i] = derivative(dx ⋅ dx)(0.0)
        end
    end

    return [Gx; Gp]
end

#-

"""
    Jacobian(f!::Function, x::Vector{Float64}, p::Float64)

Devuelve el jacobiano `J::Matrix{FLoat64}` del sistema de ecuaciones diferenciales
asociado a `f!` evaluado en `x` y `p`.
"""
function Jacobian(f!::Function, x::Vector{Float64}, p::Float64)
    s = Taylor1([0.0,0.0],1) 
    t = Taylor1(1)
    J = zeros(length(x),length(x) + 1)
    dx = [s for i in 1:length(x)]
    for i in 1:length(x)
        for j in 1:length(x)
            f!(dx,x + [j == k ? t : s for k in 1:length(x)],p)
            J[i,j] = derivative(dx[i] .+ s)(0.0)
        end
        f!(dx,x,p + t)
        J[i,end] = derivative(dx[i] .+ s)(0.0)
    end
    return J
end

#-

"""
    Jacobian(f!::Function, x::Vector{Float64}, p::Vector{Float64}, indice::Int64)

Devuelve el jacobiano `J::Matrix{FLoat64}` del sistema de ecuaciones diferenciales
asociado a `f!` evaluado en `x` y `p[indice]`.
"""
function Jacobian(f!::Function, x::Vector{Float64}, p::Vector{Float64}, indice::Int64)
    s = Taylor1([0.0,0.0],1) 
    t = Taylor1(1)
    T = [i == indice ? t : s for i in 1:length(p)]
    J = zeros(length(x),length(x) + 1)
    dx = [s for i in 1:length(x)]
    for i in 1:length(x)
        for j in 1:length(x)
            f!(dx,x + [j == k ? t : s for k in 1:length(x)],p)
            J[i,j] = derivative(dx[i] .+ s)(0.0)
        end
        f!(dx,x,p + T)
        J[i,end] = derivative(dx[i] .+ s)(0.0)
    end
    return J
end