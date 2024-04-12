# # Bifurcaciones
# Este archivo contiene las funciones correspondientes para encontrar
# puntos de bifurcación en un sistema de ecuaciones diferenciales.

#-

include("diff_tools.jl")

#-

using TaylorSeries, LinearAlgebra

#-

"""
    Bifurcation_point(f::Function,x0::Float64,p0::Float64,t::Float64; type_bif = false)

Esta función devuelve el punto de bifurcación `p::Float64, x::Float64` más cercano a `x0` y `p0`
de la ecuación diferencial descrita por `f`. Si `type_bif = true` se imprimirá el tipo de bifurcación.
"""
function Bifurcation_point(f::Function,x0::Float64,p0::Float64,t::Float64; type_bif = false)
    S = set_variables("s",numvars = 2,order = 3)
    x = x0 + S[1]
    p = p0 + S[2]
    ceros = zeros(length(x0)+1)
    Bif_vec(x,p,t) = [f(x,p,t),derivative(f(x,p,t),1)]
    i = 1
    while i <= 30 && norm(Bif_vec(x,p,t)(ceros)) > 1.e-16
        x, p = [x,p] - inv(TaylorSeries.jacobian(Bif_vec(x,p,t),ceros))*Bif_vec(x,p,t)(ceros)
        i += 1
    end
    if type_bif
        Type_Bifurcation(f,x,p,t,ceros)
    end
    return p(ceros), x(ceros)
end

#-

"""
    Bifurcation_point(f::Function,x0::Float64,p0::Float64,t::Float64,indice::Int64; type_bif = false)

Esta función devuelve el punto de bifurcación `p[indice]::Float64, x::Float64` más cercano a `x0` y `p0[indice]`
de la ecuación diferencial descrita por `f`. Si `type_bif = true` se imprimirá el tipo de bifurcación.
"""
function Bifurcation_point(f::Function,x0::Float64,p0::Vector{Float64},t::Float64,indice::Int64; type_bif = false)
    S = set_variables("s",numvars = 2,order = 3)
    x = x0 + S[1]
    p = p0 + [i == indice ? S[2] : p0[i] for i in 1:length(p)]
    p_i = p0[indice]
    ceros = zeros(length(x0)+1)
    Bif_vec(x,p,t) = [f(x,p,t),derivative(f(x,p,t),1)]
    i = 1
    while i <= 30 && norm(Bif_vec(x,p,t)(ceros)) > 1.e-16
        x, p_i = [x,p_i] - inv(TaylorSeries.jacobian(Bif_vec(x,p,t),ceros))*Bif_vec(x,p,t)(ceros)
        p = p + [i == indice ? p_i : p0[i] for i in 1:length(p)]
        i += 1
    end
    if type_bif
        Type_Bifurcation(f,x,p,t,ceros)
    end
    return p(ceros)[indice], x(ceros)
end

#-

"""
    Hopf_bifurcation(f!::Function,t::Float64,equilibrium_points::Tuple{Vector{Float64}, Vector{Vector{Float64}}})

Esta función devuelve un arreglo de los puntos `hopf_bifs::Vector{Union{Vector{Float64},Float64}}` de la rama de equilibrio 
`equilibrium_points` más cercanos a alguna bifurcación de Hopf.
"""
function Hopf_bifurcation(f!::Function,t::Float64,equilibrium_points::Tuple{Vector{Float64}, Vector{Vector{Float64}}})

    X = equilibrium_points[2]
    P = equilibrium_points[1]

    J = Jacobian.(f!,X,P,t)
    J = [j[:,1:end-1] for j in J]

    Eigvalues = [real.(eigvals(j)) for j in J]

    hopf_bifs = Vector{Union{Vector{Float64},Float64}}[]

    for i in 1:(length(Eigvalues)-1)
        j = 0
        for k in 1:length(Eigvalues[i])
            if real(Eigvalues[i][k]) * real(Eigvalues[i+1][k]) <= 0.0
                j += 1
            end
        end
    
        if j >= 2
            push!(hopf_bifs,[P[i],X[i]])
        end
    end

    if length(hopf_bifs) == 0
        return push!(hopf_bifs,[NaN,[NaN for i in X[1]]])
    end

    return hopf_bifs
end

#-

"""
    Hopf_bifurcation(f!::Function,p_ini::Vector{Float64},t::Float64,indice::Int64,
    equilibrium_points::Tuple{Vector{Float64}, Vector{Vector{Float64}}})

Esta función devuelve un arreglo de los puntos `hopf_bifs::Vector{Union{Vector{Float64},Float64}}` de la rama de equilibrio 
`equilibrium_points` más cercanos a alguna bifurcación de Hopf.
"""
function Hopf_bifurcation(f!::Function,p_ini::Vector{Float64},t::Float64,indice::Int64,equilibrium_points::Tuple{Vector{Float64}, Vector{Vector{Float64}}})

    X = equilibrium_points[2]
    P = equilibrium_points[1]
    Ps = [[i == indice ? P[j] : p_ini[i] for i in 1:length(p_ini)] for j in 1:length(P)]

    J = Jacobian.(f!,X,Ps,t,indice)
    J = [j[:,1:end-1] for j in J]

    Eigvalues = [real.(eigvals(j)) for j in J]

    hopf_bifs = Vector{Union{Vector{Float64},Float64}}[]

    for i in 1:(length(Eigvalues)-1)
        j = 0
        for k in 1:length(Eigvalues[i])
            if real(Eigvalues[i][k]) * real(Eigvalues[i+1][k]) <= 0.0
                j += 1
            end
        end
    
        if j >= 2
            push!(hopf_bifs,[P[i],X[i]])
        end
    end

    return hopf_bifs
end

#-

"""
    Type_Bifurcation(f::Function,x::TaylorN,p::TaylorN,t::Float64,ceros::Vector{Float64})

Esta es una función auxiliar para determinar el tipo de bifurcación.
"""
function Type_Bifurcation(f::Function,x::TaylorN,p::TaylorN,t::Float64,ceros::Vector{Float64})
    f_p = derivative(f(x,p,t),2)(ceros)
    f_xx = derivative(derivative(f(x,p,t),1),1)(ceros)
    f_pp = derivative(derivative(f(x,p,t),2),2)(ceros)
    f_xp = derivative(derivative(f(x,p,t),1),2)(ceros)
    f_xpp = derivative(derivative(derivative(f(x,p,t),1),2),2)(ceros)
    f_xxx = derivative(derivative(derivative(f(x,p,t),1),1),1)(ceros)

    if abs(f_p) > 1.e-16  && abs(f_xx) > 1.e-16
        println("Esta es una bifurcación Nodo de Saddle")
    elseif abs(f_p) <= 1.e-16  && abs(f_xx) > 1.e-16
        println("Esta es una bifurcación Transcrítica")
    elseif abs(f_p) <= 1.e-16  && abs(f_pp) <= 1.e-16 && abs(f_xp) > 1.e-16 && abs(f_xpp) > 1.e-16 && abs(f_xxx) > 1.e-16
        println("Esta es una bifurcación de Pitchfork")
    else
        println("No se encontró el tipo de bifurcación, estas son las derivadas parciales en ese punto:")
        println("∂f/∂p = $(f_p)")
        println("∂²f/∂x² = $(f_xx)")
        println("∂²f/∂p² = $(f_pp)")
        println("∂³f/∂x∂p² = $(f_xpp)")
        println("∂³f/∂x³= $(f_xxx)")
    end
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