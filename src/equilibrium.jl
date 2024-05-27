# # Puntos de Equilibrio
#
# En este archivo se contienen las funciones que generan las ramas de
# equilibrio de un sistema de ecuaciones diferenciales.

#-

#using TaylorSeries, LinearAlgebra

#-

#include("equilibrium_functions.jl")

#-

# #### Caso de variable unidimensinal con parametros escalares o vectoriales, permitiendo en este último escoger el parámetro de la familia de soluciones

#-

"""
    Equilibrium(f::Function, x_ini::Float64, p_ini::Float64, Δs::Float64, p_fin::Float64; N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Float64})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 
`f(x(p),p) ≈ 0.0`.
"""
function Equilibrium(f::Function, x_ini::Float64, p_ini::Float64, t::Float64, Δs::Float64, p_fin::Float64; N = 10000)
        
    X = [x_ini]
    P = [p_ini]

    x = x_ini
    p = p_ini

    i = 1

    x0 = x
    p0 = p
    x0_s, p0_s = Equilibrium_initial_values(f,x0,p0,t,p_fin)

    while p_ini <= p <= p_fin && i <= N
        
        j = 1
        while j <= 30 && norm(Equilibrium_Function(f,x,p,t,x0,p0,x0_s,p0_s,Δs)) > 1.e-16
            x, p = [x; p] - inv(Equilibrium_Jacobian(f,x,p,t,x0_s,p0_s))*Equilibrium_Function(f,x,p,t,x0,p0,x0_s,p0_s,Δs)
            j += 1
        end
        push!(X,x)
        push!(P,p)

        x0_s, p0_s = Equilibrium_Jacobian(f,x,p,t,x0_s,p0_s)\[0.0 ; 1.0]
        x0 = x
        p0 = p
        i += 1
    end

    return P, X

end

#-

"""
    Equilibrium(f::Function, x_ini::Float64, p_ini::Vector{Float64}, Δs::Float64, p_fin::Float64,indice::Int64; N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Float64})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 
`f(x(p),p') ≈ 0.0` donde `p'[indice] = p`.
"""
function Equilibrium(f::Function, x_ini::Float64, p_ini::Float64, t::Float64, Δs::Float64, p_fin::Float64,indice::Int64; N = 10000)
        
    X = [x_ini]
    P = [p_ini]

    x = x_ini
    p = p_ini

    i = 1

    x0 = x
    p0 = p
    x0_s, p0_s = Equilibrium_initial_values(f,x0,p0,t,p_fin,indice)

    while p_ini <= p <= p_fin && i <= N
        
        j = 1
        while j <= 30 && norm(Equilibrium_Function(f,x,p,t,x0,p0,x0_s,p0_s,Δs,p_ini,indice)) > 1.e-16
            x, p = [x; p] - inv(Equilibrium_Jacobian(f,x,p,t,x0_s,p0_s,p_ini,indice))*Equilibrium_Function(f,x,p,t,x0,p0,x0_s,p0_s,Δs,p_ini,indice)
            j += 1
        end
        push!(X,x)
        push!(P,p)

        x0_s, p0_s = Equilibrium_Jacobian(f,x,p,t,x0_s,p0_s,p_ini,indice)\[0.0 ; 1.0]
        x0 = x
        p0 = p
        i += 1
    end

    return P, X

end

#-

# #### Caso de variable vectorial con parametros escalares o vectoriales, permitiendo en este último escoger el parámetro de la familia de soluciones

#-

"""
    Equilibrium(f!::Function, x_ini::Vector{Float64}, p_ini::Float64, Δs::Float64, p_fin::Float64; tol = 1.e-16, N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Vector{Float64}})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 

        f!(dx,x(p),p)
        norm(dx) ≈ 0.0

"""
function Equilibrium(f!::Function, x_ini::Vector{Float64}, p_ini::Float64, t::Float64, Δs::Float64, p_fin::Float64; N = 10000)
        
    n = length(x_ini)
    X = [x_ini]
    P = [p_ini]

    x = x_ini
    p = p_ini

    i = 1

    x0 = x
    p0 = p
    x0_s, p0_s = Equilibrium_initial_values(f!,x0,p0,t,p_fin)

    while p_ini <= p <= p_fin && i <= N
        
        j = 1
        while j <= 30 && norm(Equilibrium_Function(f!,x,p,t,x0,p0,x0_s,p0_s,Δs)) > 1.e-16
            V = [x; p] - inv(Equilibrium_Jacobian(f!,x,p,t,x0_s,p0_s))*Equilibrium_Function(f!,x,p,t,x0,p0,x0_s,p0_s,Δs)
            x = V[1:n]
            p = V[end]
            j += 1
        end
        push!(X,x)
        push!(P,p)

        Dir = Equilibrium_Jacobian(f!,x,p,t,x0_s,p0_s)\[zeros(n) ; 1.0]
        x0_s = Dir[1:n]
        p0_s = Dir[end]
        x0 = x
        p0 = p
        i += 1
    end

    return P, X

end

#-

"""
    Equilibrium(f!::Function, x_ini::Vector{Float64}, p_ini::Float{Float64}, Δs::Float64, p_fin::Float64,indice::Int64; N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Vector{Float64}})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 

        f!(dx,x(p),p')
        norm(dx) ≈ 0.0

donde `p'[indice] = p`
"""
function Equilibrium(f!::Function, x_ini::Vector{Float64}, p_ini::Vector{Float64}, t::Float64, Δs::Float64, p_fin::Float64,indice::Int64; N = 10000)
        
    n = length(x_ini)
    X = [x_ini]
    P = [p_ini[indice]]

    x0_s, p0_s = Equilibrium_initial_values(f!,x_ini,p_ini,t,p_fin,indice)

    x = x_ini
    p = p_ini[indice]
    x0 = x
    p0 = p
    
    i = 1

    while p_ini[indice] <= p <= p_fin && i <= N
        
        j = 1
        while j <= 30 && norm(Equilibrium_Function(f!,x,p,t,x0,p0,x0_s,p0_s,Δs,p_ini,indice)) > 1.e-16
            V = [x; p] - inv(Equilibrium_Jacobian(f!,x,p,t,x0_s,p0_s,p_ini,indice))*Equilibrium_Function(f!,x,p,t,x0,p0,x0_s,p0_s,Δs,p_ini,indice)
            x = V[1:n]
            p = V[end]
            j += 1
        end
        push!(X,x)
        push!(P,p)

        Dir = Equilibrium_Jacobian(f!,x,p,t,x0_s,p0_s,p_ini,indice)\[zeros(n) ; 1.0]
        x0_s = Dir[1:n]
        p0_s = Dir[end]
        x0 = x
        p0 = p
        i += 1
    end

    return P, X

end

#-

# ### Referencias

#-

# - Doedel, E.J. (2007). Lecture Notes on Numerical Analysis of Nonlinear Equations.
#   In: Krauskopf, B., Osinga, H.M., Galán-Vioque, J. (eds) Numerical Continuation 
#   Methods for Dynamical Systems. Understanding Complex Systems. Springer, Dordrecht. 
#   https://doi.org/10.1007/978-1-4020-6356-5_1
