# # Funciones de puntos de Hopf
# Este archivo contiene funciones correspondientes al sistema
# de Hopf así como su jacobiano y los valores iniciales

#-

#using TaylorSeries, LinearAlgebra

#-

function Hopf_initial_values(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},t::Float64)

    n = length(equilibrium_branch[2][1])
    N_points = length(equilibrium_branch[1])
    s = Taylor1(1)
    J_old = zeros(n,n)
    Indices = Int64[]
    Eigenvals = []
    Eigenspace = []
    NO_Eigenspace = []

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

        v = [1.0,1.0]

        if !all((real.(λ_new) .* real.(λ_old)) .> 0.0) && !all((imag.(λ_new) .* imag.(λ_old)) .<= 0.0) && i > 2
            if all((imag.(λ_old) .^ 2) .< (imag.(λ_new) .^ 2))
                push!(Indices,i-1)
                push!(Eigenvals,imag(λ_old[1])^2)
                for i in 1:10
                    w = (J_old^2 + imag(λ_old[1])^2*I(n))\v
                    v = w/norm(w)
                end
                push!(Eigenspace,v)
                push!(NO_Eigenspace,nullspace(v')[:,1])

            else
                push!(Indices,i)
                push!(Eigenvals,imag(λ_new[1])^2)
                for i in 1:10
                    w = (J_new^2 + imag(λ_new[1])^2*I(n))\v
                    v = w/norm(w)
                end
                push!(Eigenspace,v)
                push!(NO_Eigenspace,nullspace(v')[:,1])
            end
        end

        J_old = J_new

    end

    return [[equilibrium_branch[1][Indices[i]],equilibrium_branch[2][Indices[i]],Eigenvals[i],Eigenspace[i],NO_Eigenspace[i]] for i in 1:length(Indices)]
end

#-

function Hopf_initial_values(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},t::Float64,p_ini::Vector{Float64},indice::Int64)

    n = length(equilibrium_branch[2][1])
    N_points = length(equilibrium_branch[1])
    s = Taylor1(1)
    r = Taylor1([0.0,0.0],3)
    J_old = zeros(n,n)
    Indices = Int64[]
    Eigenvals = []
    Eigenspace = []
    NO_Eigenspace = []

    for i in 1:n
        for j in 1:n
            dx = [s for i in 1:n]
            f!(dx,equilibrium_branch[2][1]+[j == k ? s : r for k in 1:n],[k == indice ? equilibrium_branch[1][1] : p_ini[k] for k in 1:length(p_ini)],t)
            J_old[i,j] = derivative(dx[j])(0.0)
        end
    end
    
    for i in 2:(N_points)
        J_new = zeros(n,n)
        for j in 1:n
            for k in 1:n
                dx = [s for i in 1:n]
                f!(dx,equilibrium_branch[2][i]+[k == l ? s : r for l in 1:n],[k == indice ? equilibrium_branch[1][i] : p_ini[k] for k in 1:length(p_ini)],t)
                J_new[j,k] = derivative(dx[j])(0.0)
            end
        end

        λ_new = eigvals(J_new)
        λ_old = eigvals(J_old)

        v = ones(n)

        if !all((real.(λ_new) .* real.(λ_old)) .> 0.0) && !all((imag.(λ_new) .* imag.(λ_old)) .<= 0.0) && i > 2
            if all((imag.(λ_old) .^ 2) .< (imag.(λ_new) .^ 2))
                push!(Indices,i-1)
                push!(Eigenvals,imag(λ_old[1])^2)
                for i in 1:10
                    w = (J_old^2 + imag(λ_old[1])^2*I(n))\v
                    v = w/norm(w)
                end
                push!(Eigenspace,v)
                push!(NO_Eigenspace,nullspace(v')[:,1])

            else
                push!(Indices,i)
                push!(Eigenvals,imag(λ_new[1])^2)
                for i in 1:10
                    w = (J_new^2 + imag(λ_new[1])^2*I(n))\v
                    v = w/norm(w)
                end
                push!(Eigenspace,v)
                push!(NO_Eigenspace,nullspace(v')[:,1])
            end
        end

        J_old = J_new

    end

    return [[equilibrium_branch[1][Indices[i]],equilibrium_branch[2][Indices[i]],Eigenvals[i],Eigenspace[i],NO_Eigenspace[i]] for i in 1:length(Indices)]
end

#-

function Hopf_Function(f!::Function,x0::Vector{Float64},p0::Float64,v0::Vector{Float64},κ0::Float64,t::Float64,w0::Vector{Float64})

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = p0 + S[end]

    dx = S[1:n]

    f!(dx,x,p,t)

    F1 = dx(ceros)
    F2 = [sum([(derivative(dx[i],l)*derivative(dx[l],k))(ceros)*v0[k] for l in 1:n, k in 1:n]) + κ0*v0[i] for i in 1:n]
    F3 = transpose(v0)*v0 - 1.0
    F4 = transpose(w0)*v0
    return [F1; F2; F3; F4]
end

#-

function Hopf_Function(f!::Function,x0::Vector{Float64},p0::Float64,v0::Vector{Float64},κ0::Float64,t::Float64,w0::Vector{Float64},p_ini::Vector{Float64},indice::Int64)

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = [i == indice ? (p0+S[end]) : p_ini[i] for i in 1:length(p_ini)]

    dx = S[1:n]

    f!(dx,x,p,t)

    F1 = dx(ceros)
    F2 = [sum([(derivative(dx[i],l)*derivative(dx[l],k))(ceros)*v0[k] for l in 1:n, k in 1:n]) + κ0*v0[i] for i in 1:n]
    F3 = transpose(v0)*v0 - 1.0
    F4 = transpose(w0)*v0
    return [F1; F2; F3; F4]
end

#-

function Hopf_Jacobian(f!::Function,x0::Vector{Float64},p0::Float64,v0::Vector{Float64},κ0::Float64,t::Float64,w0::Vector{Float64})

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = p0 + S[end]

    dx = S[1:n]

    f!(dx,x,p,t)

    A = [derivative(dx[i],j)(ceros) for i in 1:n, j in 1:n]
    B = [derivative(dx[i],n+1)(ceros) for i in 1:n]
    C = [0.0 for i in 1:n, j in 1:n]
    D = zeros(n)
    
    E = [sum([derivative(derivative(dx[i],l)*derivative(dx[l],k)*v0[k],j)(ceros) for  l in 1:n, k in 1:n]) for i in 1:n, j in 1:n]
    F = [sum([derivative(derivative(dx[i],l)*derivative(dx[l],k)*v0[k],n+1)(ceros) for  l in 1:n, k in 1:n]) for i in 1:n]
    G = [sum([(derivative(dx[i],k)*derivative(dx[k],j))(ceros) for k in 1:n]) + ==(i,j)*κ0 for i in 1:n, j in 1:n]
    K = v0

    L = zeros(1,n)
    M = 0.0
    N = 2*transpose(v0)
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

function Hopf_Jacobian(f!::Function,x0::Vector{Float64},p0::Float64,v0::Vector{Float64},κ0::Float64,t::Float64,w0::Vector{Float64},p_ini::Vector{Float64},indice::Int64)

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = [i == indice ? (p0+S[end]) : p_ini[i] for i in 1:length(p_ini)]

    dx = S[1:n]

    f!(dx,x,p,t)

    A = [derivative(dx[i],j)(ceros) for i in 1:n, j in 1:n]
    B = [derivative(dx[i],n+1)(ceros) for i in 1:n]
    C = [0.0 for i in 1:n, j in 1:n]
    D = zeros(n)
    
    E = [sum([derivative(derivative(dx[i],l)*derivative(dx[l],k)*v0[k],j)(ceros) for  l in 1:n, k in 1:n]) for i in 1:n, j in 1:n]
    F = [sum([derivative(derivative(dx[i],l)*derivative(dx[l],k)*v0[k],n+1)(ceros) for  l in 1:n, k in 1:n]) for i in 1:n]
    G = [sum([(derivative(dx[i],k)*derivative(dx[k],j))(ceros) for k in 1:n]) + ==(i,j)*κ0 for i in 1:n, j in 1:n]
    K = v0

    L = zeros(1,n)
    M = 0.0
    N = 2*transpose(v0)
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

