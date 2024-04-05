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
function Gradient(f::Function,x::Float64,p::Float64,t::Float64)
    return [derivative(f(x+s,p+r,t+r))(0.0),derivative(f(x+r,p+s,t+r))(0.0)]
end

#-

"""
    Gradient(f::Function,x::Float64,p::Vector{Float64},indice::Int64)

Devuelve el gradiente `∇f::Vector{Float64}` de la función `f` evaluado
en `x` y `p[indice]`.
"""
function Gradient(f::Function,x::Float64,p::Vector{Float64},t::Float64,indice::Int64)
    return [derivative(f(x+s,p .+ r,t+r))(0.0); [i == indice ? derivative(f(x+r,p+S,t+r))(0.0) : 0.0 for i in 1:length(p)]]
end

#-

"""
    Gradient(f::Function,x::Vector{Float64},p::Float64)

Devuelve el gradiente `∇||f||::Vector{Float64}` de la norma del sistema
de ecuaciones diferenciales asociado a función `f!` evaluado en `x` y `p`.
"""
function Gradient(f!::Function,x::Vector{Float64},p::Float64,t::Float64)
    dx = [r for i in 1:length(x)]
    G = zeros(length(x) + 1)
    for i in 1:length(x)
        f!(dx,x + [i == j ? s : r for j in 1:length(x)],p + r, t + r)
        G[i] = derivative(dx ⋅ dx)(0.0)
    end
    f!(dx,x .+ r,p + s, t + r)
    G[end] = derivative(dx ⋅ dx)(0.0)
    return G
end

#-

"""
    Gradient(f!::Function,x::Vector{Float64},p::Vector{Float64},indice::Int64)

Devuelve el gradiente `∇||f||::Vector{Float64}` de la norma del sistema
de ecuaciones diferenciales asociado a función `f!` evaluado en `x` y `p[indice]`.
"""
function Gradient(f!::Function,x::Vector{Float64},p::Vector{Float64},t::Float64,indice::Int64)

    S = [i == indice ? s : r for i in 1:length(p_ini)]

    dx = [r for i in 1:length(x)]
    Gx = zeros(length(x))
    Gp = zeros(length(p))

    for i in 1:length(x)
        f!(dx,x + [i == j ? s : r for j in 1:length(x)],p .+ r, t + r)
        Gx[i] = derivative(dx ⋅ dx)(0.0)
    end

    for i in 1:length(p)
        f!(dx,x .+ r,p + S, t + r)
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
function Jacobian(f!::Function, x::Vector{Float64}, p::Float64,t::Float64)

    J = zeros(length(x),length(x) + 1)
    dx = [r for i in 1:length(x)]
    for i in 1:length(x)
        for j in 1:length(x)
            f!(dx,x + [j == k ? s : r for k in 1:length(x)],p + r, t + r)
            J[i,j] = derivative(dx[i])(0.0)
        end
        f!(dx,x .+ r,p + s,t + r)
        J[i,end] = derivative(dx[i])(0.0)
    end
    return J
end

#-

"""
    Jacobian(f!::Function, x::Vector{Float64}, p::Vector{Float64}, indice::Int64)

Devuelve el jacobiano `J::Matrix{FLoat64}` del sistema de ecuaciones diferenciales
asociado a `f!` evaluado en `x` y `p[indice]`.
"""
function Jacobian(f!::Function, x::Vector{Float64}, p::Vector{Float64},t::Float64, indice::Int64)
    S = [i == indice ? s : r for i in 1:length(p_ini)]
    J = zeros(length(x),length(x) + 1)
    dx = [s for i in 1:length(x)]
    for i in 1:length(x)
        for j in 1:length(x)
            f!(dx,x + [j == k ? s : r for k in 1:length(x)],p .+ r,t + r)
            J[i,j] = derivative(dx[i])(0.0)
        end
        f!(dx,x .+ r,p + S,t + r)
        J[i,end] = derivative(dx[i] )(0.0)
    end
    return J
end

function Hessian(f!::Function, x::Vector{Float64}, p::Float64,t::Float64)
    S = set_variables("s", numvars = (length(x) + 1), order = 2)
    dx = x .+ S[1:end-1]
    H = zeros(length(x) + 1,length(x) + 1)
    for i in 1:length(x) + 1
        for j in 1:length(x) + 1
            f!(dx,x .+ S[1:end-1],p + S[end],t)
            H[i,j] = derivative(derivative(dx ⋅ dx,i),j)(zeros(length(x) + 1))
        end
    end
    return H
end