# # Puntos de Equilibrio
#
# En este archivo se contienen las funciones que generan las ramas de
# equilibrio de un sistema de ecuaciones diferenciales.

#-

using TaylorSeries, LinearAlgebra

#-

s = Taylor1(1)
r = Taylor1([0.0,0.0],1)

#-

include("diff_tools.jl")
include("implicit_function.jl")
include("newton.jl")
include("psarc.jl")

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
        
    X = Float64[]
    P = Float64[]

    push!(P,p_ini)
    push!(X,x_ini)

    x_s, p_s = Derivative_arclength(f, x_ini, p_ini, t, p_fin)
    x = x_ini + x_s*Δs
    p = p_ini + p_s*Δs

    x, p = Newton(f,x,p,t)

    push!(P,p)
    push!(X,x)

    i = 0

    while p_ini <= p <= p_fin && i <= N 
        x_s, p_s = step(f, x+x_s*Δs, p+p_s*Δs, t, x_s, p_s)
        x, p = Newton(f,x + Δs*x_s,p + Δs*p_s, t)
        push!(P,p)
        push!(X,x)
        i+=1  
    end

    return P,X

end

#-

"""
    Equilibrium(f::Function, x_ini::Float64, p_ini::Vector{Float64}, Δs::Float64, p_fin::Float64,indice::Int64; N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Float64})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 
`f(x(p),p') ≈ 0.0` donde `p'[indice] = p`.
"""
function Equilibrium(f::Function, x_ini::Float64, p_ini::Vector{Float64}, t::Float64, Δs::Float64, p_fin::Float64,indice::Int64; N = 10000)
    
    S = [i == indice ? s : r for i in 1:length(p_ini)]

    X = Float64[]
    P = Float64[]

    if indice > length(p_ini)
        throw(ArgumentError("El índice de variable no está dentro de la dimensión"))
    end

    push!(P,p_ini[indice])
    push!(X,x_ini)

    x_s, p_s = Derivative_arclength(f, x_ini, p_ini, t, p_fin, indice)
    x = x_ini + x_s*Δs
    p = p_ini + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p_ini)]

    x, p = Newton(f,x,p,t,indice)
    push!(P,p[indice])
    push!(X,x)
    

    i = 0

    while p_ini[indice] <= p[indice] <= p_fin && i <= N
        x_s, p_s = step(f!, x + x_s*Δs, p + p_s*Δs, t, x_s, p_s,indice)
        x, p = Newton(f,x + Δs*x_s,p + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p_ini)], t,indice)
        push!(P,p[indice])
        push!(X,x)
        i+=1  
    end

    return P,X

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
        
    X = Vector{Float64}[]
    P = Float64[]

    push!(P,p_ini)
    push!(X,x_ini)

    x_s, p_s = Derivative_arclength(f!, x_ini, p_ini, t, p_fin)
    x_s, p_s = step(f!, x_ini + x_s*Δs, p_ini + p_s*Δs, t, x_s, p_s)

    x = x_ini + x_s*Δs
    p = p_ini + p_s*Δs

    
    x, p = Newton(f!,x,p,t)

    push!(P,p)
    push!(X,x)
    
    i = 2

    while p_ini <= p <= p_fin && i <= N 
        x_s, p_s = step(f!, x + x_s*Δs, p + p_s*Δs, t, x_s, p_s)
        x, p = Newton(f!,x + x_s*Δs, p + p_s*Δs, t)
        push!(P,p)
        push!(X,x)
        i+=1  
    end

    return P,X

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
function Equilibrium(f!::Function, x_ini::Vector{Float64}, p_ini::Vector{Float64}, t::Float64, Δs::Float64, p_fin::Float64,indice::Int64; tol = 1.e-16, N = 10000)
    
    S = [i == indice ? s : r for i in 1:length(p_ini)]

    X = Vector{Float64}[]
    P = Float64[]

    if indice > length(p_ini)
        throw(ArgumentError("El índice de variable no está dentro de la dimensión"))
    end

    push!(P,p_ini[indice])
    push!(X,x_ini)

    x_s, p_s = Derivative_arclength(f!, x_ini, p_ini, t, p_fin,indice)

    x = x_ini + x_s*Δs
    p = p_ini + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p_ini)]
    
    x, p = Newton(f!,x,p,t,indice)

    push!(P,p[indice])
    push!(X,x)

    i = 0

    while p_ini[indice] <= p[indice] <= p_fin && i <= N
        x_s, p_s = step(f!, x + x_s*Δs, p + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p_ini)], t, x_s, p_s,indice)
        x, p = Newton(f!,x + x_s*Δs,p + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p_ini)], t,indice)
        push!(P,p[indice])
        push!(X,x)
        i+=1  
    end

    return P,X

end

#-

# ### Referencias

#-

# - Doedel, E.J. (2007). Lecture Notes on Numerical Analysis of Nonlinear Equations.
#   In: Krauskopf, B., Osinga, H.M., Galán-Vioque, J. (eds) Numerical Continuation 
#   Methods for Dynamical Systems. Understanding Complex Systems. Springer, Dordrecht. 
#   https://doi.org/10.1007/978-1-4020-6356-5_1
