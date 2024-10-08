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
function Equilibrium(f::Function, x_ini::Float64, p_ini::Float64, t::Float64, Δs::Float64, N::Int64; newtonite = 8, newtontol = 1.e-16, bifurcations = false)
    
    #

    x = zeros(N)
    p = zeros(N)
    stable = [false for i in 1:N]
    unstable = [true for i in 1:N]

    if bifurcations
        Fold_Index = Int64[]
    end

    s = Taylor1(1)
    r = zero(s)

    F = zeros(2)
    J = zeros(2,2)

    xs = 0.0
    ps = 0.0

    Δ = zeros(2)
    ΔS = zeros(2)

    dx_old = 0.0

    #

    x[1] = x_ini
    p[1] = p_ini

    Equilibrium_Jacobian!(J, f, x[1], p[1], t, xs, ps, s, r)

    ΔS .= nullspace([J[1, 1] J[1, 2]])
    xs = ΔS[1]
    ps = ΔS[2]

    if J[1,1] <= 0
        stable[1] = true
        unstable[1] = false
    end

    dx_old = J[1,1]

    #

    Equilibrium_Function!(F, f, x[1], p[1], t, x[1], p[1], xs, ps, Δs)
    Equilibrium_Jacobian!(J, f, x[1], p[1], t, xs, ps, s, r)

    for i in 2:N

        x[i] = x[i-1]
        p[i] = p[i-1]
        
        j = 1

        while j <= newtonite && norm(F) > newtontol

            Δ .= -inv(J)*F
            x[i] += Δ[1]
            p[i] += Δ[2]

            Equilibrium_Function!(F, f, x[i], p[i], t, x[i-1], p[i-1], xs, ps, Δs)
            Equilibrium_Jacobian!(J, f, x[i], p[i], t, xs, ps, s, r)

            j += 1

        end

        ΔS .= [J[1,1] J[1,2]; xs ps] \ [0.0; 1.0]

        xs = ΔS[1]/norm(ΔS)
        ps = ΔS[2]/norm(ΔS)

        Equilibrium_Function!(F, f, x[i], p[i], t, x[i-1], p[i-1], xs, ps, Δs)
        Equilibrium_Jacobian!(J, f, x[i], p[i], t, xs, ps, s, r)

        if J[1,1] <= 0
            stable[i] = true
            unstable[i] = false
        end

        if bifurcations
    
            if dx_old*J[1,1] < 0
                push!(Fold_Index, i)
            end

            dx_old = J[1,1]
    
        end

    end

    if bifurcations

        pf = p[Fold_Index]
        xf = x[Fold_Index]

        Limit_Points!(f, pf, xf; newtonite = newtonite, newtontol = newtontol)

        return p, x, stable, unstable, pf, xf

    else
        return p, x, stable, unstable
    end

end

#-

"""
    Equilibrium(f::Function, x_ini::Float64, p_ini::Float64, Δs::Float64, p_fin::Float64; N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Float64})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 
`f(x(p),p) ≈ 0.0`.
"""
function Equilibrium(f::Function, x_ini::Float64, p_ini::Vector{Float64}, t::Float64, Δs::Float64, indice::Int64, N::Int64; newtonite = 8, newtontol = 1.e-16)
    
    #

    np = length(p_ini)

    x = zeros(N)
    p = zeros(N, np)
    stable = [false for i in 1:N]
    unstable = [true for i in 1:N]

    s = Taylor1(1)
    r = zero(s)

    Sp = [i == indice ? s : r for i in 1:np]

    F = zeros(2)
    J = zeros(2,2)

    xs = 0.0
    ps = 0.0

    Δ = zeros(2)
    ΔS = zeros(2)

    #

    x[1] = x_ini
    p[1, :] .= p_ini

    Equilibrium_Jacobian!(J, f, x[1], p[1, :], t, xs, ps, s, r, Sp)

    ΔS .= nullspace([J[1, 1] J[1, 2]])
    xs = ΔS[1]
    ps = ΔS[2]

    if J[1,1] <= 0
        stable[1] = true
        unstable[1] = false
    end

    #

    Equilibrium_Function!(F, f, x[1], p[1, :], t, x[1], p[1, :], xs, ps, Δs, indice)
    Equilibrium_Jacobian!(J, f, x[1], p[1, :], t, xs, ps, s, r, Sp)

    for i in 2:N

        x[i] = x[i-1]
        p[i, :] .= p[i-1, :]
        
        j = 1

        while j <= newtonite && norm(F) > newtontol

            Δ .= -inv(J)*F

            x[i] += Δ[1]
            p[i, indice] += Δ[2]

            Equilibrium_Function!(F, f, x[i], p[i, :], t, x[i-1], p[i-1, :], xs, ps, Δs, indice)
            Equilibrium_Jacobian!(J, f, x[i], p[i, :], t, xs, ps, s, r, Sp)

            j += 1

        end

        ΔS .= [J[1,1] J[1,2]; xs ps] \ [0.0; 1.0]

        xs = ΔS[1]
        ps = ΔS[2]

        Equilibrium_Function!(F, f, x[i], p[i, :], t, x[i-1], p[i-1, :], xs, ps, Δs, indice)
        Equilibrium_Jacobian!(J, f, x[i], p[i, :], t, xs, ps, s, r, Sp)

        if J[1,1] <= 0
            stable[i] = true
            unstable[i] = false
        end
        
    end

    return p, x, stable, unstable

end

#-

"""
    Equilibrium(f::Function, x_ini::Float64, p_ini::Float64, Δs::Float64, p_fin::Float64; N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Float64})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 
`f(x(p),p) ≈ 0.0`.
"""
function Equilibrium(f!::Function, x_ini::Vector{Float64}, p_ini::Float64, t::Float64, Δs::Float64, N::Int64; newtonite = 8, newtontol = 1.e-16, bifurcations = false)
    
    #

    n = length(x_ini)

    x = zeros(N, n)
    p = zeros(N)

    stable = [false for i in 1:N]
    unstable = [true for i in 1:N]

    if bifurcations
        Fold_Index = Int64[]
        Hopf_Index = Int64[]
        Fold_Test = Int64[]
        Hopf_Test = Int64[]
    end
    
    max_norm = 0.0
    max_ind = 0

    λ = zeros(Complex{Float64}, n)
    λ_old = zeros(Complex{Float64}, n)

    s = Taylor1(1)
    r = zero(s)

    dx = zeros(n)
    dxT = [s for i in 1:n]

    F = zeros(n + 1)
    J = zeros(n + 1,n + 1)

    Sx = [i == j ? s : r for i in 1:n, j in 1:n]

    xs = zeros(n)
    ps = 0.0

    Δ  = zeros(n + 1)
    ΔS = zeros(n + 1)

    #

    x[1, :] .= x_ini
    p[1]     = p_ini

    Equilibrium_Jacobian!(J, f!, x[1, :], p[1], t, xs, ps, dxT, s, r, Sx, n)

    ΔS .= nullspace(J[1:n, :])
    xs .= ΔS[1:n]
    ps  = ΔS[n+1]

    λ .= eigvals(J[1:n, 1:n])
    λ_old .= λ

    hopf = [real(λ[i]) * real(λ_old[i]) < 0.0 && imag(λ[i]) * imag(λ_old[i]) > 0.0 for i in 1:n]
    fold = [real(λ[i]) * real(λ_old[i]) < 0.0 && imag(λ[i]) * imag(λ_old[i]) < 0.0 for i in 1:n]

    if all(real.(λ) .< 0.0)
        stable[1] = true
        unstable[1] = false
    end

    #

    Equilibrium_Function!(F, f!, x[1, :], p[1], t, x[1, :], p[1], xs, ps, Δs, dx, n)
    Equilibrium_Jacobian!(J, f!, x[1, :], p[1], t, xs, ps, dxT, s, r, Sx, n)

    normF = norm(F[1:n])

    if normF > max_norm
        max_norm = normF
        max_ind = 1
    end

    for i in 2:N

        x[i, :] .= x[i-1, :] .+ Δs*xs
        p[i]     = p[i-1]     + Δs*ps
        
        j = 1

        while j <= newtonite && norm(F) > newtontol

            Δ .= -inv(J)*F
            x[i, :] .+= Δ[1:n]
            p[i] += Δ[end]

            Equilibrium_Function!(F, f!, x[i, :], p[i], t, x[i-1, :], p[i-1], xs, ps, Δs, dx, n)
            Equilibrium_Jacobian!(J, f!, x[i, :], p[i], t, xs, ps, dxT, s, r, Sx, n)

            j += 1

        end

        normF = norm(F[1:n])

        if normF > max_norm
            max_norm = normF
            max_ind = 1
        end

        ΔS .= J \ [zeros(n); 1.0]

        xs .= ΔS[1:n]/norm(ΔS)
        ps = ΔS[end]/norm(ΔS)

        Equilibrium_Function!(F, f!, x[i, :], p[i], t, x[i-1, :], p[i-1], xs, ps, Δs, dx, n)
        Equilibrium_Jacobian!(J, f!, x[i, :], p[i], t, xs, ps, dxT, s, r, Sx, n)

        λ .= eigvals(J[1:n, 1:n])

        if all(real.(λ) .< 0.0)
            stable[i] = true
            unstable[i] = false
        end
    
        if bifurcations

            hopf .= [real(λ[i]) * real(λ_old[i]) <= 0.0 && imag(λ[i]) * imag(λ_old[i]) > 0.0 for i in 1:n]
            fold .= [real(λ[i]) * real(λ_old[i]) <= 0.0 && imag(λ[i]) * imag(λ_old[i]) <= 0.0 for i in 1:n]
    
            if any(hopf)
                push!(Hopf_Index, i)
                push!(Hopf_Test, (1:n)[hopf][1])
            end

            if any(fold)
                push!(Fold_Index, i)
                push!(Fold_Test, (1:n)[fold][1])
            end
    
            λ_old .= λ
    
        end

    end

    if bifurcations

        pf = [p[i] for i in Fold_Index]
        xf = [x[i,j] for i in Fold_Index, j in 1:n]

        pb = [p[i] for i in Hopf_Index]
        xb = [x[i,j] for i in Hopf_Index, j in 1:n]

        Limit_Points!(f!, pf, xf, Fold_Test; newtonite = newtonite, newtontol = newtontol)
        Hopf_Points!(f!, pb, xb, Hopf_Test; newtonite = newtonite, newtontol = newtontol)

        println("┌ Datos de las rams de equilibrio")

        println("├ Número de puntos: $(N)")

        println("├ Norma máxima del sistema: $(max_norm) (punto # $(max_ind))")

        for i in 1:length(Fold_Index)
            Equilibrium_Jacobian!(J, f!, xf[i, :], pf[i], t, xs, ps, dxT, s, r, Sx, n)
            λ = eigvals(J[1:n,1:n])
            println("├ Bifurcacion LP $(i):")
            for i in 1:n
                println("├ \t λ_$(i) = $(λ[i])")
            end
        end

        for i in 1:length(Hopf_Index)
            Equilibrium_Jacobian!(J, f!, xb[i, :], pb[i], t, xs, ps, dxT, s, r, Sx, n)
            λ = eigvals(J[1:n,1:n])
            println("├ Bifurcacion Hopf $(i):")
            for i in 1:n
                println("├ \t λ_$(i) = $(λ[i])")
            end
        end

        println("└")

        return p, x, stable, unstable, pf, xf, pb, xb
    else

        println("┌ Datos de las rams de equilibrio")

        println("├ Número de puntos: $(N)")

        println("└ Norma máxima del sistema: $(max_norm) (punto # $(max_ind))")

        return p, x, stable, unstable
    end

end

#-

"""
    Equilibrium(f::Function, x_ini::Float64, p_ini::Float64, Δs::Float64, p_fin::Float64; N = 10000)

Devuelve una tupla de vectores `(P::Vector{Float64}, X::Vector{Float64})` donde `P` 
es el dominio de la rama de solución correspondiente a las condiciones iniciales 
dadas y `X` es la imagen los valores de `P`, entendiendo que `x(p)` con `x ∈ X` y `p ∈ P` es tal que 
`f(x(p),p) ≈ 0.0`.
"""
function Equilibrium(f!::Function, x_ini::Vector{Float64}, p_ini::Vector{Float64}, t::Float64, Δs::Float64, indice::Int64, N::Int64; newtonite = 8, newtontol = 1.e-16, bifurcations = false)
    
    #

    n = length(x_ini)
    np = length(p_ini)

    x = zeros(N, n)
    p = [p_ini[j] for i in 1:N, j in 1:np]

    stable = [false for i in 1:N]
    unstable = [true for i in 1:N]

    if bifurcations
        Fold_Index = Int64[]
        Hopf_Index = Int64[]
        Fold_Test = Int64[]
        Hopf_Test = Int64[]
    end

    max_norm = 0.0
    max_ind = 0

    normF = 0.0

    λ = zeros(Complex{Float64}, n)
    λ_old = zeros(Complex{Float64}, n)

    s = Taylor1(1)
    r = zero(s)

    dx = zeros(n)
    dxT = [s for i in 1:n]

    F = zeros(n + 1)
    J = zeros(n + 1,n + 1)

    Sx = [i == j ? s : r for i in 1:n, j in 1:n]
    Sp = [i == indice ? s : r for i in 1:np]

    xs = zeros(n)
    ps = 0.0

    Δ  = zeros(n + 1)
    ΔS = zeros(n + 1)

    V = [zeros(n); 1.0]

    #

    x[1, :] .= x_ini
    #p[1, :] .= p_ini

    Equilibrium_Jacobian!(J, f!, x[1, :], p[1, :], t, xs, ps, dxT, s, r, Sx, Sp, n)

    ΔS .= nullspace(J[1:n, :])

    for j in 1:n 
        xs[j] = ΔS[j]
    end

    ps = ΔS[n+1]

    λ .= eigvals(J[1:n, 1:n])
    λ_old .= λ

    hopf = [real(λ[i]) * real(λ_old[i]) <= 0.0 && imag(λ[i]) * imag(λ_old[i]) > 0.0 for i in 1:n]
    fold = [real(λ[i]) * real(λ_old[i]) <= 0.0 && imag(λ[i]) * imag(λ_old[i]) <= 0.0 for i in 1:n]

    if all(real.(λ) .< 0.0)
        stable[1] = true
        unstable[1] = false
    end

    #

    Equilibrium_Function!(F, f!, x[1, :], p[1, :], t, x[1, :], p[1, :], xs, ps, Δs, dx, n, indice)
    Equilibrium_Jacobian!(J, f!, x[1, :], p[1, :], t, xs, ps, dxT, s, r, Sx, Sp, n)

    normF = norm(F[1:n])

    if normF > max_norm
        max_norm = normF
        max_ind = 1
    end

    for i in 2:N
        
        for j in 1:n 
            x[i, j] = x[i-1, j] + Δs*xs[j]
        end

        p[i, indice] = p[i-1, indice] + Δs*ps
        
        j = 1

        while j <= newtonite && norm(F) > newtontol

            Δ .= J \ F

            for k in 1:n
                x[i, k] -= Δ[k]
            end

            p[i, indice] -= Δ[n+1]

            Equilibrium_Function!(F, f!, x[i, :], p[i, :], t, x[i-1, :], p[i-1, :], xs, ps, Δs, dx, n, indice)
            Equilibrium_Jacobian!(J, f!, x[i, :], p[i, :], t, xs, ps, dxT, s, r, Sx, Sp, n)

            j += 1

        end

        normF = norm(F[1:n])

        if normF > max_norm
            max_norm = normF
            max_ind = i
        end

        ΔS .= J \ V

        for j in 1:n
            xs[j] = ΔS[j]/norm(ΔS)
        end

        ps  = ΔS[n+1]/norm(ΔS)

        Equilibrium_Function!(F, f!, x[i, :], p[i, :], t, x[i-1, :], p[i-1, :], xs, ps, Δs, dx, n, indice)
        Equilibrium_Jacobian!(J, f!, x[i, :], p[i, :], t, xs, ps, dxT, s, r, Sx, Sp, n)

        λ .= eigvals(J[1:n, 1:n])

        if all(real.(λ) .< 0.0)
            stable[i] = true
            unstable[i] = false
        end
    
        if bifurcations

            hopf .= [real(λ[i]) * real(λ_old[i]) <= 0.0 && imag(λ[i]) * imag(λ_old[i]) > 0.0 for i in 1:n]
            fold .= [real(λ[i]) * real(λ_old[i]) <= 0.0 && imag(λ[i]) * imag(λ_old[i]) <= 0.0 for i in 1:n]
    
            if any(hopf)
                push!(Hopf_Index, i)
                push!(Hopf_Test, (1:n)[hopf][1])
            end

            if any(fold)
                push!(Fold_Index, i)
                push!(Fold_Test, (1:n)[fold][1])
            end
    
            λ_old .= λ
    
        end

    end

    if bifurcations

        pf = [p[i,j] for i in Fold_Index, j in 1:np]
        xf = [x[i,j] for i in Fold_Index, j in 1:n]

        pb = [p[i,j] for i in Hopf_Index, j in 1:np]
        xb = [x[i,j] for i in Hopf_Index, j in 1:n]

        Limit_Points!(f!, pf, xf, Fold_Test, indice; newtonite = newtonite, newtontol = newtontol)
        Hopf_Points!(f!, pb, xb, Hopf_Test, indice; newtonite = newtonite, newtontol = newtontol)

        println("┌ Datos de las rams de equilibrio")

        println("├ Número de puntos: $(N)")

        println("├ Norma máxima del sistema: $(max_norm) (punto # $(max_ind))")

        for i in 1:length(Fold_Index)
            Equilibrium_Jacobian!(J, f!, xf[i, :], pf[i, :], t, xs, ps, dxT, s, r, Sx, Sp, n)
            λ = eigvals(J[1:n,1:n])
            println("├ Bifurcacion LP $(i):")
            for i in 1:n
                println("├ \t λ_$(i) = $(λ[i])")
            end
        end

        for i in 1:length(Hopf_Index)
            Equilibrium_Jacobian!(J, f!, xb[i, :], pb[i, :], t, xs, ps, dxT, s, r, Sx, Sp, n)
            λ = eigvals(J[1:n,1:n])
            println("├ Bifurcacion Hopf $(i):")
            for i in 1:n
                println("├ \t λ_$(i) = $(λ[i])")
            end
        end

        println("└")

        return p, x, stable, unstable, pf, xf, pb, xb
    else

        println("┌ Datos de las rams de equilibrio")

        println("├ Número de puntos: $(N)")

        println("└ Norma máxima del sistema: $(max_norm) (punto # $(max_ind))")

        return p, x, stable, unstable
    end
end

