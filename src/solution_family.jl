# # Familias de Soluciones
#
# En este archivo se contienen las funciones que generan las familias de
# soluciones de un sistema unidimensional o multidimensional.

using TaylorSeries, LinearAlgebra

include("psarc.jl")
include("newton.jl")

# # Caso de variable unidimensinal con parametros escalares o vectoriales, permitiendo en este último escoger el parámetro de la familia de soluciones

"""
    Solution_family(f::Function, x_ini::Float64, p_ini::Float64, Δs::Float64, p_fin::Float64,orden::Int64; tol = 1.e-16, N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Float64})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 
`f(x(p),p) ≈ 0.0`.
"""
function Solution_family(f::Function, x_ini::Float64, p_ini::Float64, Δs::Float64, p_fin::Float64,orden::Int64; tol = 1.e-16, N = 10000)
        
    X = Float64[]
    P = Float64[]

    t = Taylor1(orden)

    push!(P,p_ini)
    push!(X,x_ini)

    x_s, p_s = first_step(f, x_ini, p_ini, p_fin,orden)
    x = x_ini + x_s*Δs
    p = p_ini + p_s*Δs

    
    x = Newton(f,x + t,p)
    p = Newton(f,x,p + t)
    if abs(f(x,p)) <= 1.e-16
        push!(P,p)
        push!(X,x)
    else
        throw(ArgumentError("No se pudo calcular el primer paso, prueba con valores iniciales distintos"))
    end
    

    i = 0

    while p_ini <= p <= p_fin && i <= N
        x_s, p_s = step(f, x, p, x_s, p_s ,orden)
        x = Newton(f,x + x_s*Δs + t,p + p_s*Δs)
        p = Newton(f,x,p + p_s*Δs + t)
        push!(P,p)
        push!(X,x)
        i+=1  
    end

    return P[1:end-1],X[1:end-1]

end

"""
    Solution_family(f::Function, x_ini::Float64, p_ini::Vector{Float64}, Δs::Float64, p_fin::Float64,orden::Int64,indice::Int64; tol = 1.e-16, N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Float64})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 
`f(x(p),p') ≈ 0.0` donde `p'[indice] = p`.
"""
function Solution_family(f::Function, x_ini::Float64, p_ini::Vector{Float64}, Δs::Float64, p_fin::Float64,orden::Int64,indice::Int64; tol = 1.e-16, N = 10000)
        
    X = Float64[]
    P = Float64[]

    if indice > length(p_ini)
        throw(ArgumentError("El índice de variable no está dentro de la dimensión"))
    end

    t = Taylor1(orden)
    T = [i == indice ? Taylor1(orden) : Taylor1(0) for i in 1:length(p_ini)]

    push!(P,p_ini[indice])
    push!(X,x_ini)

    x_s, p_s = first_step(f, x_ini, p_ini, p_fin,orden,indice)
    x = x_ini + x_s*Δs
    p = p_ini + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p_ini)]

    
    x = Newton(f,x + t,p)
    p = Newton(f,x,p + T,indice)
    if abs(f(x,p)) <= 1.e-16
        push!(P,p[indice])
        push!(X,x)
    else
        throw(ArgumentError("No se pudo calcular el primer paso, prueba con valores iniciales distintos"))
    end
    

    i = 0

    while p_ini[indice] <= p[indice] <= p_fin && i <= N
        x_s, p_s = step(f, x, p, x_s, p_s ,orden,indice)
        x = Newton(f,x + x_s*Δs + t,p + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p)])
        p = Newton(f,x,p + [i == indice ? p_s*Δs + t : Taylor1(0) for i in 1:length(p)],indice)
        push!(P,p[indice])
        push!(X,x)
        i+=1  
    end

    return P[1:end-1],X[1:end-1]

end

# # Caso de variable vectorial con parametros escalares o vectoriales, permitiendo en este último escoger el parámetro de la familia de soluciones

"""
    Solution_family(f!::Function, x_ini::Vector{Float64}, p_ini::Float64, Δs::Float64, p_fin::Float64,orden::Int64; tol = 1.e-16, N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Vector{Float64}})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 

        f!(dx,x(p),p)
        norm(dx) ≈ 0.0

"""
function Solution_family(f!::Function, x_ini::Vector{Float64}, p_ini::Float64, Δs::Float64, p_fin::Float64,orden::Int64; tol = 1.e-16, N = 10000)
        
    X = Vector{Float64}[]
    P = Float64[]

    t = Taylor1(orden)

    push!(P,p_ini)
    push!(X,x_ini)

    x_s, p_s = first_step(f!, x_ini, p_ini, p_fin,orden)

    #println("p_s = $(p_s)")
    
    x = x_ini + x_s*Δs
    p = p_ini + p_s*Δs

    
    x = Newton(f!,x,p,orden)
    p = Newton(f!,x,p + t)

    dx_old = [Taylor1(0) for i in 1:length(x_ini)]
    f!(dx_old,x,p)

    #println("dx = $(dx_old(0.0))")
    #println("|dx| = $(norm(dx_old(0.0)))")

    if norm(dx_old(0.0)) <= 1.e-16
        push!(P,p)
        push!(X,x)
    else
        throw(ArgumentError("No se pudo calcular el primer paso, prueba con valores iniciales distintos"))
    end
    

    i = 0

    while p_ini <= p <= p_fin && i <= N

        #println("x_s = $(x_s)")
        #println("x = $(x)")

        x_s, p_s = step(f!, x, p, x_s, p_s ,orden)

        #println("x_s = $(x_s)")
        
        x = Newton(f!,x + x_s*Δs,p + p_s*Δs,orden)
        p = Newton(f!,x,p + p_s*Δs + t)
        push!(P,p)
        push!(X,x)
        i+=1  
    end

    return P[1:end-1],X[1:end-1]

end


"""
    Solution_family(f!::Function, x_ini::Vector{Float64}, p_ini::Float{Float64}, Δs::Float64, p_fin::Float64,orden::Int64,indice::Int64; tol = 1.e-16, N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Vector{Float64}})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 

        f!(dx,x(p),p')
        norm(dx) ≈ 0.0

donde `p'[indice] = p`
"""
function Solution_family(f!::Function, x_ini::Vector{Float64}, p_ini::Vector{Float64}, Δs::Float64, p_fin::Float64,orden::Int64,indice::Int64; tol = 1.e-16, N = 10000)
        
    X = Vector{Float64}[]
    P = Float64[]

    if indice > length(p_ini)
        throw(ArgumentError("El índice de variable no está dentro de la dimensión"))
    end

    t = Taylor1(orden)
    T = [i == indice ? Taylor1(orden) : Taylor1(0) for i in 1:length(p_ini)]

    push!(P,p_ini[indice])
    push!(X,x_ini)

    x_s, p_s = first_step(f!, x_ini, p_ini, p_fin,orden,indice)

    x = x_ini + x_s*Δs
    p = p_ini + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p_ini)]
    
    x = Newton(f!,x,p,orden)
    p = Newton(f!,x,p + T,orden,indice)
    
    dx_old = [Taylor1(0) for i in 1:length(x_ini)]
    f!(dx_old,x,p)

    if norm(dx_old(0.0)) <= 1.e-16
        push!(P,p[indice])
        push!(X,x)
    else
        throw(ArgumentError("No se pudo calcular el primer paso, prueba con valores iniciales distintos"))
    end
    

    i = 0

    while p_ini[indice] <= p[indice] <= p_fin && i <= N
        x_s, p_s = step(f!, x, p, x_s, p_s ,orden,indice)
        x = Newton(f!,x + x_s*Δs,p + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p)],orden)
        p = Newton(f!,x,p + [i == indice ? p_s*Δs + t : Taylor1(0) for i in 1:length(p)],orden,indice)
        push!(P,p[indice])
        push!(X,x)
        i+=1  
    end

    return P[1:end-1],X[1:end-1]

end

