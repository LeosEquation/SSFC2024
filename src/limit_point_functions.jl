# # Funciones de puntos límite
# Este archivo contiene funciones correspondientes al sistema
# de punto límite así como su jacobiano y los valores iniciales

#-

#using TaylorSeries, LinearAlgebra

#-

function LP_initial_values(f::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Float64}},t::Float64)

    N_points = length(equilibrium_branch[1])
    s = Taylor1(1)
    J_old = 0.0
    Indices = Int64[]
    Eigenvecs = Float64[]
    Eigenvals = Float64[]

    J_old = derivative(f(equilibrium_branch[2][1]+s,equilibrium_branch[1][1],t))(0.0)
    
    for i in 2:(N_points)

        J_new = derivative(f(equilibrium_branch[2][i]+s,equilibrium_branch[1][i],t))(0.0)

        λ_new = J_new
        λ_old = J_old

        if (λ_new * λ_old) < 0.0 && i > 2
            v_old = eigvecs(J_old)[:,1][1]
            v_new = eigvecs(J_new)[:,1][1]
            if abs(λ_old) < abs(λ_new)
                push!(Indices,i-1)
                push!(Eigenvals,λ_old)
                push!(Eigenvecs,v_old)
            else
                push!(Indices,i)
                push!(Eigenvals,λ_new)
                push!(Eigenvecs,v_new)
            end
        end

        J_old = J_new

    end

    return [[equilibrium_branch[1][Indices[i]],equilibrium_branch[2][Indices[i]],Eigenvecs[i],Eigenvals[i]] for i in 1:length(Indices)]
end

#-

function LP_initial_values(f::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Float64}},t::Float64,p_ini::Vector{Float64},indice::Int64)

    N_points = length(equilibrium_branch[1])
    s = Taylor1(1)
    J_old = 0.0
    Indices = Int64[]
    Eigenvecs = Float64[]
    Eigenvals = Float64[]

    for i in 1:n
        for j in 1:n
            J_old = derivative(f(equilibrium_branch[2][1]+s,[k == indice ? equilibrium_branch[1][1] : p_ini[k] for k in 1:length(p_ini)],t))(0.0)
        end
    end
    
    for i in 2:(N_points)
        J_new = 0.0
        for j in 1:n
            for k in 1:n
                J_new = derivative(derivative(f(equilibrium_branch[2][i]+s,[k == indice ? equilibrium_branch[1][i] : p_ini[k] for k in 1:length(p_ini)],t)))(0.0)
            end
        end

        λ_new = J_new
        λ_old = J_old

        if (λ_new * λ_old) > 0.0 && i > 2
            v_old = eigvecs(J_old)[:,1]
            v_new = eigvecs(J_new)[:,1]
            if abs(λ_old) < abs(λ_new)
                push!(Indices,i-1)
                push!(Eigenvals,λ_old)
                push!(Eigenvecs,v_old)
            else
                push!(Indices,i)
                push!(Eigenvals,λ_new)
                push!(Eigenvecs,v_new)
            end
        end

        J_old = J_new

    end

    return [[equilibrium_branch[1][Indices[i]],equilibrium_branch[2][Indices[i]],Eigenvecs[i],Eigenvals[i]] for i in 1:length(Indices)]
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

        if !all((real.(λ_new) .* real.(λ_old)) .> 0.0) && !all((imag.(λ_new) .* imag.(λ_old)) .> 0.0) && i > 2
            norm_new = findmin(norm.(λ_new))
            norm_old = findmin(norm.(λ_old))
            v_old = eigvecs(J_old)[:,norm_old[2]]
            v_new = eigvecs(J_new)[:,norm_new[2]]
            if norm_old[1] < norm_new[1] && all(imag(v_old) .== 0.0) && imag(norm_old[1]) == 0.0
                push!(Indices,i-1)
                push!(Eigenvals,norm_old[1])
                push!(Eigenvecs,v_old)
            elseif norm_old[1] >= norm_new[1] && all(imag(v_new) .== 0.0) && imag(norm_new[1]) == 0.0
                push!(Indices,i)
                push!(Eigenvals,norm_new[1])
                push!(Eigenvecs,v_new)
            end
        end

        J_old = J_new

    end

    return [[equilibrium_branch[1][Indices[i]],equilibrium_branch[2][Indices[i]],Eigenvecs[i],Eigenvals[i]] for i in 1:length(Indices)]
end

#-

function LP_initial_values(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},t::Float64,p_ini::Vector{Float64},indice::Int64)

    n = length(equilibrium_branch[2][1])
    N_points = length(equilibrium_branch[1])
    s = Taylor1(1)
    r = Taylor1([0.0,0.0],3)
    J_old = zeros(n,n)
    Indices = Int64[]
    Eigenvecs = Vector[]
    Eigenvals = []

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

        if !all((real.(λ_new) .* real.(λ_old)) .> 0.0) && !all((imag.(λ_new) .* imag.(λ_old)) .> 0.0) && i > 2
            norm_new = findmin(norm.(λ_new))
            norm_old = findmin(norm.(λ_old))
            v_old = eigvecs(J_old)[:,norm_old[2]]
            v_new = eigvecs(J_new)[:,norm_new[2]]
            if norm_old[1] < norm_new[1] && all(imag(v_old) .== 0.0) && imag(norm_old[1]) == 0.0
                push!(Indices,i-1)
                push!(Eigenvals,norm_old[1])
                push!(Eigenvecs,v_old)
            elseif norm_old[1] >= norm_new[1] && all(imag(v_new) .== 0.0) && imag(norm_new[1]) == 0.0
                push!(Indices,i)
                push!(Eigenvals,norm_new[1])
                push!(Eigenvecs,v_new)
            end
        end

        J_old = J_new

    end

    return [[equilibrium_branch[1][Indices[i]],equilibrium_branch[2][Indices[i]],real(Eigenvecs[i]),real(Eigenvals[i])] for i in 1:length(Indices)]
end

#-

function LP_Jacobian(f::Function,x0::Float64,p0::Float64,v0::Float64,t::Float64; ω = 1.0)

    ceros = zeros(2)

    S = set_variables("s", numvars = 2, order = 3)

    x = x0 + S[1]
    p = p0 + S[2]

    dx = S[1]

    dx = f(x,p,t)
    
    A = derivative(dx,1)(ceros)
    B = derivative(dx,2)(ceros)
    C = 0.0

    D = derivative(derivative(dx,1),1)(ceros)*v0
    E = derivative(derivative(dx,1),2)(ceros)*v0
    F = A
    
    G = 0.0
    K = 0.0
    M = ω
    
    Extended_Jacobian = [A B C;
                         D E F;
                         G K M]
    
    return Extended_Jacobian
end

#-

function LP_Jacobian(f::Function,x0::Float64,p0::Float64,v0::Float64,t::Float64,p_ini::Vector{Float64},indice::Int64; ω = 1.0)

    ceros = zeros(2)

    S = set_variables("s", numvars = 2, order = 3)

    x = x0 + S[1]
    p = [i == indice ? (p0 + S[2]) : p_ini[i] for i in 1:length(p_ini)]

    dx = S[1]

    dx = f(x,p,t)
    
    A = derivative(dx,1)(ceros)
    B = derivative(dx,2)(ceros)
    C = 0.0

    D = derivative(derivative(dx,1),1)(ceros)*v0
    E = derivative(derivative(dx,1),2)(ceros)*v0
    F = A
    
    G = 0.0
    K = 0.0
    M = ω
    
    Extended_Jacobian = [A B C;
                         D E F;
                         G K M]
    
    return Extended_Jacobian
end

#-

function LP_Jacobian(f!::Function,x0::Vector{Float64},p0::Float64,v0::Union{Vector{Float64},Vector{ComplexF64}},t::Float64; ω = ones(length(x0)))

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

    D = [sum([derivative(derivative(dx[i],k)j)(ceros)*v0[k] for k in 1:n]) for i in 1:n, j in 1:n]
    E = [sum([derivative(derivative(dx[i],k)n+1)(ceros)*v0[k] for k in 1:n]) for i in 1:n]
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

function LP_Jacobian(f!::Function,x0::Vector{Float64},p0::Float64,v0::Union{Vector{Float64},Vector{ComplexF64}},t::Float64, p_ini::Vector{Float64}, indice::Int64; ω = ones(length(x0)))

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

    D = [sum([derivative(derivative(dx[i],k)j)(ceros)*v0[k] for k in 1:n]) for i in 1:n, j in 1:n]
    E = [sum([derivative(derivative(dx[i],k)n+1)(ceros)*v0[k] for k in 1:n]) for i in 1:n]
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

function LP_Function(f::Function,x0::Float64,p0::Float64,v0::Float64,t::Float64; ω = 1.0)

    ceros = zeros(2)

    S = set_variables("s", numvars = 2, order = 3)

    x = x0 + S[1]
    p = p0 + S[2]

    dx = f(x,p,t)

    LP1 = dx(ceros)
    LP2 = derivative(dx,1)(ceros)*v0
    LP3 = ω*v0 - 1.0

    return [LP1; LP2; LP3]
end

#-

function LP_Function(f::Function,x0::Float64,p0::Float64,v0::Float64,t::Float64,p_ini::Vector{Float64}, indice::Int64; ω = ones(length(x0)))

    ceros = zeros(2)

    S = set_variables("s", numvars = 2, order = 3)

    x = x0 + S[1]
    p = [i == indice ? (p0 + S[2]) : p_ini[i] for i in 1:length(p_ini)]

    dx = S[1]

    dx = f(x,p,t)

    LP1 = dx(ceros)
    LP2 = derivative(dx,1)(ceros)*v0
    LP3 = ω*v0 - 1.0

    return [LP1; LP2; LP3]
end

#-

function LP_Function(f!::Function,x0::Vector{Float64},p0::Float64,v0::Union{Vector{Float64},Vector{ComplexF64}},t::Float64; ω = ones(length(x0)))

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = p0 + S[end]

    dx = S[1:n]

    f!(dx,x,p,t)

    LP1 = dx(ceros)
    LP2 = [sum([derivative(dx[i],k)(ceros)*v0[k] for k in 1:n]) for i in 1:n]
    LP3 = ω ⋅ v0 - 1.0

    return [LP1; LP2; LP3]
end

#-

function LP_Function(f!::Function,x0::Vector{Float64},p0::Float64,v0::Union{Vector{Float64},Vector{ComplexF64}},t::Float64, p_ini::Vector{Float64}, indice::Int64; ω = ones(length(x0)))

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = [i == indice ? (p0+S[end]) : p_ini[i] for i in 1:length(p_ini)]

    dx = S[1:n]

    f!(dx,x,p,t)

    LP1 = dx(ceros)
    LP2 = [sum([derivative(dx[i],k)(ceros)*v0[k] for k in 1:n]) for i in 1:n]
    LP3 = ω ⋅ v0 - 1.0

    return [LP1; LP2; LP3]
end