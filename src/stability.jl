# # Estabilidad
#
# En este archivo se contienen las funciones que determinan intervalos
# y zonas de estabilidad en las ramas de equilibrio del sistema de 
# ecuaciones diferenciales.

#-

#using TaylorSeries

#-

"""
    Stability(f::Function,x::Float64,p::Float64,t::Float64)

Esta función devuelve `true` si el punto de equilibrio situado en `p::Float64, x::Float64`
es estable y `false` en caso contrario.
"""
function Stability(f::Function,x::Float64,p::Float64,t::Float64)
    if derivative(f(x+Taylor1(1),p,t)) <= 0.0
        return true
    else
        return false
    end
end

#-

"""
    Stability_intervals(f::Function,t::Float64,equilibrium_points::Tuple{Vector{Float64}, Vector{Float64}})

Esta función devuelve 2 arreglos `stability_intervals::Vector{Tuple}, unstability_intervals::Vector{Tuple}` donde 
`stability_intervals` consta de tuplas `(p::Float64,x::Float64)` con valores `p` y `x` siendo puntos estables 
sobre la rama de soluciones `equilibrium_points` y tuplas `(NaN,NaN)` cuando los valores sobre la rama son inestables.
COntrariamente, `stability_intervals` consta de tuplas `(x,p)` cuando son puntos inestables sobre `equilibrium_points`
y `(NaN,NaN)` en caso contrario. 

La salida de esta función está enfocada a facilitar la graficación de zonas de estabilidad e inestabilidad simultaneamente.
"""
function Stability_intervals(f::Function,t::Float64,equilibrium_points::Tuple{Vector{Float64}, Vector{Float64}})
    
    X = equilibrium_points[2]
    P = equilibrium_points[1]
    stability_intervals = [(P[i],X[i]) for i in 1:length(P)]
    unstability_intervals = [(P[i],X[i]) for i in 1:length(P)]
    
    for i in 1:length(P)
        if Stability(f,X[i],P[i],t)
            unstability_intervals[i] = (NaN,NaN)
        else
            stability_intervals[i] = (NaN,NaN)
        end
    end
    return stability_intervals, unstability_intervals
end

#-

"""
    Stability(f::Function,x::Float64,p::Vector{Float64},t::Float64)

Esta función devuelve `true` si el punto de equilibrio situado en `p::Vector{Float64}, x::Float64`
es estable y `false` en caso contrario.
"""
function Stability(f::Function,x::Float64,p::Vector{Float64},t::Float64)
    if derivative(f(x+Taylor1(1),p,t)) <= 0.0
        return true
    else
        return false
    end
end

#-

"""
    Stability_intervals(f::Function,p_ini::Vector{Float64},t::Float64,indice::Int64,
    equilibrium_points::Tuple{Vector{Float64}, Vector{Float64}})

Esta función devuelve 2 arreglos `stability_intervals::Vector{Tuple}, unstability_intervals::Vector{Tuple}` donde 
`stability_intervals` consta de tuplas `(p::Float64,x::Float64)` con valores `p` y `x` siendo puntos estables 
sobre la rama de soluciones `equilibrium_points` y tuplas `(NaN,NaN)` cuando los valores sobre la rama son inestables.
COntrariamente, `stability_intervals` consta de tuplas `(x,p)` cuando son puntos inestables sobre `equilibrium_points`
y `(NaN,NaN)` en caso contrario. 

La salida de esta función está enfocada a facilitar la graficación de zonas de estabilidad e inestabilidad simultaneamente.
"""
function Stability_intervals(f::Function,p_ini::Vector{Float64},t::Float64,indice::Int64,equilibrium_points::Tuple{Vector{Float64}, Vector{Float64}})
    
    X = equilibrium_points[2]
    P = equilibrium_points[1]
    Ps = [[j == indice ? P[i] : p_ini[j] for j in 1:length(p_ini)] for i in 1:length(P)]

    stability_intervals = [(P[i],X[i]) for i in 1:length(P)]
    unstability_intervals = [(P[i],X[i]) for i in 1:length(P)]
    
    for i in 1:length(P)
        if Stability(f,X[i],Ps[i],t)
            unstability_intervals[i] = (NaN,NaN)
        else
            stability_intervals[i] = (NaN,NaN)
        end
    end
    return stability_intervals, unstability_intervals
end

#-

"""
    Stability(f::Function,x::Vector{Float64},p::Float64,t::Float64)

Esta función devuelve `true` si el punto de equilibrio situado en `p::Float64, x::Vector{Float64}`
es estable y `false` en caso contrario.
"""
function Stability(f!::Function,x::Vector{Float64},p::Float64,t::Float64)

    n = length(x)
    s = Taylor1(2)
    r = Taylor1([0.0,0.0],2)

    J = zeros(n,n)
    for i in 1:n
        for j in 1:n
            dx = [s for i in 1:n]
            f!(dx,x + [k == j ? s : r for k in 1:n], p + r, t + r)
            J[i,j] = differentiate(dx[i])(0.0)
        end
    end

    if all(real.(eigvals(J)) .< 0.0)
        return true
    else
        return false
    end
end

#-

"""
    Stability_intervals(f::Function,t::Float64,equilibrium_points::Tuple{Vector{Float64}, Vector{Vector{Float64}}})

Esta función devuelve 2 arreglos `stability_intervals::Vector{Tuple}, unstability_intervals::Vector{Tuple}` donde 
`stability_intervals` consta de tuplas `(p::Float64,x::Vector{Float64})` con valores `p` y `x` siendo puntos estables 
sobre la rama de soluciones `equilibrium_points` y tuplas `(NaN,[NaN,...])` cuando los valores sobre la rama son inestables.
Contrariamente, `stability_intervals` consta de tuplas `(x,p)` cuando son puntos inestables sobre `equilibrium_points`
y `(NaN,[NaN,...])` en caso contrario. 

La salida de esta función está enfocada a facilitar la graficación de zonas de estabilidad e inestabilidad simultaneamente.
"""
function Stability_intervals(f!::Function,t::Float64,equilibrium_points::Tuple{Vector{Float64}, Vector{Vector{Float64}}})
    
    X = equilibrium_points[2]
    P = equilibrium_points[1]

    stability_intervals = [(P[i],X[i]) for i in 1:length(P)]
    unstability_intervals = [(P[i],X[i]) for i in 1:length(P)]
    
    for i in 1:length(P)
        if Stability(f!,X[i],P[i],t)
            unstability_intervals[i] = (NaN,[NaN for j in 1:length(X[i])])
        else
            stability_intervals[i] = (NaN,[NaN for j in 1:length(X[i])])
        end
    end
    return stability_intervals, unstability_intervals
end

#-

"""
    Stability(f::Function,x::Vector{Float64},p::Vector{Float64},t::Float64)

Esta función devuelve `true` si el punto de equilibrio situado en `p::Vector{Float64}, x::Vector{Float64}`
es estable y `false` en caso contrario.
"""
function Stability(f!::Function,x::Vector{Float64},p::Vector{Float64},t::Float64,indice::Int64)

    n = length(x)
    s = Taylor1(2)
    r = Taylor1([0.0,0.0],2)

    J = zeros(n,n)
    for i in 1:n
        for j in 1:n
            dx = [s for i in 1:n]
            f!(dx,x + [k == j ? s : r for k in 1:n], p .+ r, t + r)
            J[i,j] = differentiate(dx[i])(0.0)
        end
    end

    if all(real.(eigvals(J)) .< 0.0)
        return true
    else
        return false
    end
end

#-

"""
    Stability_intervals(f!::Function,p_ini::Vector{Float64},t::Float64,indice::Int64,
    equilibrium_points::Tuple{Vector{Float64}, Vector{Vector{Float64}}})

Esta función devuelve 2 arreglos `stability_intervals::Vector{Tuple}, unstability_intervals::Vector{Tuple}` donde 
`stability_intervals` consta de tuplas `(p::Float64,x::Vector{Float64})` con valores `p` y `x` siendo puntos estables 
sobre la rama de soluciones `equilibrium_points` y tuplas `(NaN,[NaN,...])` cuando los valores sobre la rama son inestables.
Contrariamente, `stability_intervals` consta de tuplas `(x,p)` cuando son puntos inestables sobre `equilibrium_points`
y `(NaN,[NaN,...])` en caso contrario. 

La salida de esta función está enfocada a facilitar la graficación de zonas de estabilidad e inestabilidad simultaneamente.
"""
function Stability_intervals(f!::Function,p_ini::Vector{Float64},t::Float64,indice::Int64,equilibrium_points::Tuple{Vector{Float64}, Vector{Vector{Float64}}})
    
    X = equilibrium_points[2]
    P = equilibrium_points[1]
    Ps = [[i == indice ? P[j] : p_ini[i] for i in 1:length(p_ini)] for j in 1:length(P)]

    stability_intervals = [(P[i],X[i]) for i in 1:length(P)]
    unstability_intervals = [(P[i],X[i]) for i in 1:length(P)]
    
    for i in 1:length(P)
        if Stability(f!,X[i],Ps[i],t,indice)
            unstability_intervals[i] = (NaN,[NaN for j in 1:length(X[i])])
        else
            stability_intervals[i] = (NaN,[NaN for j in 1:length(X[i])])
        end
    end
    return stability_intervals, unstability_intervals
end

#-

# ### Referencias

#-

# - http://www.maths.liv.ac.uk/~bnvasiev/Past%20students/Caitlin_399.pdf
# 
# - Doedel, E.J. (2007). Lecture Notes on Numerical Analysis of Nonlinear Equations.
#   In: Krauskopf, B., Osinga, H.M., Galán-Vioque, J. (eds) Numerical Continuation 
#   Methods for Dynamical Systems. Understanding Complex Systems. Springer, Dordrecht. 
#   https://doi.org/10.1007/978-1-4020-6356-5_1 