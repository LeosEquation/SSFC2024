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

function LP_initial_values(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},t::Float64)

    n = length(equilibrium_branch[2][1])
    N_points = length(equilibrium_branch[1])
    s = Taylor1(1)
    J_old = zeros(n,n)
    Indices = Int64[]
    Eigenvecs = Vector[]
    Eigenvals = []

    for i in 1:n
        for j in 1:n
            dx = [s for i in 1:n]
            f!(dx,equilibrium_branch[2][1]+[j == k ? s : 0.0 for k in 1:n],equilibrium_branch[1][1],t)
            J_old[i,j] = derivative(dx[j])(0.0)
        end
    end
    
    for i in 2:(N_points)
        J_new = zeros(n,n)
        for j in 1:n
            for k in 1:n
                dx = [s for i in 1:n]
                f!(dx,equilibrium_branch[2][i]+[k == l ? s : 0.0 for l in 1:n],equilibrium_branch[1][i],t)
                J_new[j,k] = derivative(dx[j])(0.0)
            end
        end

        λ_new = eigvals(J_new)
        λ_old = eigvals(J_old)

        if !all((real.(λ_new) .* real.(λ_old)) .> 0.0) && !all((imag.(λ_new) .* imag.(λ_old)) .> 0.0)
            norm_new = findmin(norm.(λ_new))
            norm_old = findmin(norm.(λ_old))
            if norm_old[1] < norm_new[1]
                push!(Indices,i-1)
                push!(Eigenvals,norm_old[1])
                push!(Eigenvecs,eigvecs(J_old)[:,norm_old[2]])
            else
                push!(Indices,i)
                push!(Eigenvals,norm_new[1])
                push!(Eigenvecs,eigvecs(J_new)[:,norm_new[2]])
            end
        end

        J_old = J_new

    end

    return [[equilibrium_branch[1][Indices[i]],equilibrium_branch[2][Indices[i]],Eigenvecs[i],Eigenvals[i]] for i in 1:length(Indices)]
end

#-

function LP_Jacobian(f!::Function,x0::Vector{Float64},p0::Float64,ν0::Union{Vector{Float64},Vector{ComplexF64}},t::Float64; ω = [1.0,1.0])

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = p0 + S[end]

    dx = x

    f!(dx,x,p,t)
    
    A = [derivative(dx[i],j)(ceros) for i in 1:n, j in 1:n]
    B = [derivative(dx[i],n+1)(ceros) for i in 1:n]
    C = [0.0 for i in 1:n, j in 1:n]

    D = [sum([derivative(derivative(dx[i],k)j)(ceros)*ν0[k] for k in 1:n]) for i in 1:n, j in 1:n]
    E = [sum([derivative(derivative(dx[i],k)n+1)(ceros)*ν0[k] for k in 1:n]) for i in 1:n]
    F = A
    
    G = zeros(1,n)
    K = 0.0
    M = transpose(ω)
    
    Extended_Jacobian = [A B C;
                         D E F;
                         G K M]
    
    return Extended_Jacobian
end

#-

function LP_Function(f!::Function,x0::Vector{Float64},p0::Float64,ν0::Union{Vector{Float64},Vector{ComplexF64}},t::Float64; ω = [1.0,1.0])

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = p0 + S[end]

    dx = x

    f!(dx,x,p,t)

    LP1 = dx(ceros)
    LP2 = [sum([derivative(dx[i],k)(ceros)*ν0[k] for k in 1:n]) for i in 1:n]
    LP3 = ω ⋅ ν0 - 1.0

    return [LP1; LP2; LP3]
end

#-

function Limit_Points(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},Δs::Float64,t::Float64; ω = [1.0,1.0], test = false)
    
    initial_values = LP_initial_values(f!,equilibrium_branch,t)

    LPs = Vector{Union{Float64,Vector{Float64}}}[]

    for initial_value in initial_values
        p0, x0, ν0, λ0 = initial_value

        x = x0
        p = p0
        ν = ν0

        i = 1

        while i <= 30 && norm(LP_Function(f!,x,p,ν,t; ω = ω)) > 1.e-16

            if test
                println("Paso $(i): λ = $(p)")
                println("           x = $(x)")
                println("           norm(LP_Function) = $(norm(LP_Function(f!,x,p,ν,t)))")
                println("\n")
                show(stdout, "text/plain",  LP_Function(f!,x,p,ν,t))
                println("\n")
            end

            X = [x;p;ν] - inv(LP_Jacobian(f!,x,p,ν,t; ω = ω))*LP_Function(f!,x,p,ν,t; ω = ω)
            x = X[1:length(x)]
            p = X[length(x)+1]
            ν = X[length(x)+2:end]
            i += 1
        end
        
        if norm(LP_Function(f!,x,p,ν,t; ω = ω)) < Δs
            push!(LPs,[p,x])
        end

    end

        return LPs
end

#-

function Hopf_initial_values(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},t::Float64)

    n = length(equilibrium_branch[2][1])
    N_points = length(equilibrium_branch[1])
    s = Taylor1(1)
    J_old = zeros(n,n)
    Indices = Int64[]
    Eigenvecs = []
    Eigenvals = []
    Nullspaces = []

    for i in 1:n
        for j in 1:n
            dx = [s for i in 1:n]
            f!(dx,equilibrium_branch[2][1]+[j == k ? s : 0.0 for k in 1:n],equilibrium_branch[1][1],t)
            J_old[i,j] = derivative(dx[j])(0.0)
        end
    end
    
    for i in 2:(N_points)
        J_new = zeros(n,n)
        for j in 1:n
            for k in 1:n
                dx = [s for i in 1:n]
                f!(dx,equilibrium_branch[2][i]+[k == l ? s : 0.0 for l in 1:n],equilibrium_branch[1][i],t)
                J_new[j,k] = derivative(dx[j])(0.0)
            end
        end

        λ_new = eigvals(J_new)
        λ_old = eigvals(J_old)

        if !all((real.(λ_new) .* real.(λ_old)) .> 0.0) && !all((imag.(λ_new) .* imag.(λ_old)) .<= 0.0)
            if all((imag.(λ_old) .^ 2) .< (imag.(λ_new) .^ 2))
                push!(Indices,i-1)
                push!(Eigenvecs,eigvecs(J_old))
                push!(Nullspaces,eigvecs(transpose(J_old^2 + imag(λ_old[1])^2*I(n))))
                push!(Eigenvals,λ_old)
            else
                push!(Indices,i)
                push!(Eigenvecs,eigvecs(J_new))
                push!(Nullspaces,eigvecs(transpose(J_new^2 + imag(λ_new[1])^2*I(n))))
                push!(Eigenvals,λ_new)
            end
        end

        J_old = J_new

    end

    return [[equilibrium_branch[1][Indices[i]],equilibrium_branch[2][Indices[i]],Eigenvecs[i],Eigenvals[i], Nullspaces[i]] for i in 1:length(Indices)]
end

#-

function Hopf_function(f!::Function,x0,p0,ν0,κ0,t::Float64,w0)

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = p0 + S[end]

    dx = x

    f!(dx,x,p,t)

    F1 = dx(ceros)
    F2 = [sum([(derivative(dx[i],l)*derivative(dx[l],k))(ceros)*ν0[k] for l in 1:n, k in 1:n]) + κ0*ν0[i] for i in 1:n]
    F3 = transpose(ν0)*ν0 - 1.0
    F4 = transpose(w0)*ν0
    return [F1; F2; F3; F4]
end

#-

function Hopf_Jacobian(f!,x0,p0,ν0,κ0,t,w0)

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = p0 + S[end]

    dx = x

    f!(dx,x,p,t)

    A = [derivative(dx[i],j)(ceros) for i in 1:n, j in 1:n]
    B = [derivative(dx[i],n+1)(ceros) for i in 1:n]
    C = [0.0 for i in 1:n, j in 1:n]
    D = zeros(n)
    
    E = [sum([derivative(derivative(dx[i],l)*derivative(dx[l],k)*ν0[k],j)(ceros) for  l in 1:n, k in 1:n]) for i in 1:n, j in 1:n]
    F = [sum([derivative(derivative(dx[i],l)*derivative(dx[l],k)*ν0[k],n+1)(ceros) for  l in 1:n, k in 1:n]) for i in 1:n]
    G = [sum([(derivative(dx[i],k)*derivative(dx[k],j))(ceros) for k in 1:n]) + ==(i,j)*κ0 for i in 1:n, j in 1:n]
    K = ν0

    L = zeros(1,n)
    M = 0.0
    N = 2*transpose(ν0)
    Ñ = 0.0

    O = zeros(1,n)
    P = 0.0
    Q = transpose(w0)
    R = 0.0

    Extended_Jacobian = [A B C D;
                         E F G K;
                         L M N Ñ;
                         O P Q R]
    
    return Extended_Jacobian
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