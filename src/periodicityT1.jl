#-

function Periodic_Data(N, max_norm, max_ind)
    println("┌ Datos de las ramas periódicas")

    println("├ Número de puntos: $(N)")

    println("└ Máximo valor de la integral ∫₀ᵀF(x,p)dt : $(max_norm) (punto # $(max_ind))")
end

#

function Periodic_Orbits(f!::Function, x_ini::Vector{Float64}, p_ini::Float64, t::Float64, Δs::Float64, N::Int64; integorder = 20, integtol = 1.e-20, newtonite = 8, newtontol = 1.e-16, maxtol = 1.e-8, taylororder = 1, integmaxsteps = 1000)

    ### Inicializando variables

    n  = length(x_ini)

    dx = zeros(Float64, n)

    dx0 = zeros(Float64, n)

    stable = [false for i in 1:N]

    Dfx = zeros(Float64, n,n)

    F = zeros(Float64, n+2)

    J = zeros(Float64, n+2,n+2)

    p = zeros(N)
    x = zeros(Float64, N, n)
    T = zeros(Float64, N)

    M = zeros(n,n)

    s = Taylor1(Float64, taylororder)
    r = Taylor1(zeros(Float64, taylororder + 1))

    S = set_variables("s", numvars = n+1, order = 1)

    zerosTN = zeros(Float64, n+1)

    dxT1 = [s for i in 1:n]

    SxT1 = [i == j ? s : r for i in 1:n, j in 1:n]

    Sx = S[1:n]
    Sp = S[n+1]

    Δ = zeros(Float64, n + 2)

    xs = zeros(Float64, n)
    ps = 0.0
    Ts = 0.0

    ΔS = zeros(Float64, n + 2)

    V = [zeros(Float64, n+1) ; 1.0]

    ### Cálculo de la frecuencia y periodo inicial

    JacobianT1!(Dfx, f!, x_ini, p_ini, t, dxT1, SxT1, r, n)

    ω_ini = abs(imag(eigvals(Dfx)[end]))

    T_ini = 2*π/ω_ini

    ### Cálculo de la órbita en el punto de equilibrio

    W = nullspace([-ω_ini*I(n) Dfx; Dfx ω_ini*I(n)])[:,1]

    ws = W[1:n]
    wc = W[n+1:end]

    #ϕ(t) = ws*sin(2*π*t/T_ini) + wc*cos(2*π*t/T_ini)
    dϕ(t) = ω_ini*(ws*cos(ω_ini*t) - wc*sin(ω_ini*t))

    ### Asignando puntos iniciales

    x[1,:] .= x_ini
    T[1]    = T_ini
    p[1]    = p_ini

    dx0 .= δx*dϕ(0.0)

    time, sol = taylorinteg(f!, x_ini .+ Sx, 0.0, T_ini, integorder, integtol, p_ini + Sp; maxsteps = integmaxsteps)

    for i in 1:n
        for j in 1:n 
            M[i,j] = differentiate(sol[end,i], j)(zerosTN)
        end
    end

    if all(norm.(eigvals(M)) .<= 1.0)
        stable[1] = true
    end

    if length(time) > integmaxsteps
        return p_ini, x_ini, p_ini, stable[1]
    end

    Jacobian_Periodic_Orbits!(J, f!, p_ini, T_ini, xs, ps, Ts, sol[end,:], dx, dx0, zerosTN, n)

    J[1:n, n+2] .=dϕ(0.0)

    ΔS .= nullspace(J[1:n+1, :])[:, 1]

    for j in 1:n
        xs[j] = ΔS[j]
    end
    ps = ΔS[n+1]
    Ts = ΔS[n+2]

    ### Calculando la primer órbita

    i = 2

    for j in 1:n 
          x[i,j] = x[i-1,j] + Δs*xs[j]
    end

    p[i] = p[i-1] + Δs*ps
    T[i] = T[i-1] + Δs*Ts

    time, sol = taylorinteg(f!, x[i,:] .+ Sx, 0.0, T[i], integorder, integtol, p[i] + Sp; maxsteps = integmaxsteps)

    if length(time) > integmaxsteps
        return p[1:i-1], x[1:i-1, :], T[1:i-1], stable[1]
    end

    Function_Periodic_Orbits!(F, x[i, :], p[i], T[i], x[i-1, :], p[i-1], T[i-1], xs, ps, Ts, Δs, sol[end, :](zerosTN), dx0, n)
    Jacobian_Periodic_Orbits!(J, f!, p[i], T[i], xs, ps, Ts, sol[end,:], dx, dx0, zerosTN, n)

    j = 1

    while j < newtonite && norm(F) > newtontol
        
        try 
            Δ .= J \ F
        catch
            println("Error: El Jacobiano del primer punto es no invertible, intenta cambiando el tamaño del paso y los parámetros de integración")
            return p_ini, x_ini, T_ini, stable[1]
        end

        for k in 1:n 
            x[i,k]   -= Δ[k]
        end

        p[i] -= Δ[n+1]
        T[i] -= Δ[n+2]

        time, sol = taylorinteg(f!, x[i,:] .+ Sx, 0.0, T[i], integorder, integtol, p[i] .+ Sp; maxsteps = integmaxsteps)

        if length(time) > integmaxsteps
            return p_ini, x_ini, T_ini, stable[1]
        end

        Function_Periodic_Orbits!(F, x[i, :], p[i], T[i], x[i-1, :], p[i-1], T[i-1], xs, ps, Ts, Δs, sol[end, :](zerosTN), dx0, n)
        Jacobian_Periodic_Orbits!(J, f!, p[i], T[i], xs, ps, Ts, sol[end,:], dx, dx0, zerosTN, n)

        j += 1

    end

    for i in 1:n
        for j in 1:n 
            M[i,j] = differentiate(sol[end,i], j)(zerosTN)
        end
    end

    if all(norm.(eigvals(M)) .<= 1.0)
        stable[i] = false
    end

    if norm(F) > maxtol
        return p_ini, x_ini, T_ini, stable[1]
    end

    println("$(i): ite = $(j) : norm(F) = $(norm(F))")

    # f!(dx0, x[i,:], p[i, :], 0.0)

    # Jacobian_Periodic_Orbits!(J, f!, p[i, :], T[i], xs, ps, Ts, sol[end,:], dx, dx0, zerosTN, n)

    ΔS .= J \ V

    for j in 1:n 
        xs[j] = ΔS[j]/norm(ΔS)
    end

    ps  = ΔS[n+1]/norm(ΔS)
    Ts  = ΔS[n+2]/norm(ΔS)

    f!(dx0, x[i,:], p[i], 0.0)

    ### Calculando el resto de puntos

    for i in 3:N

        for j in 1:n 
            x[i,j] = x[i-1,j] + Δs*xs[j]
        end

        p[i] = p[i-1] + Δs*ps
        T[i] = T[i-1] + Δs*Ts

        time, sol = taylorinteg(f!, x[i,:] .+ Sx, 0.0, T[i], integorder, integtol, p[i] + Sp; maxsteps = integmaxsteps)

        if length(time) > integmaxsteps
            return p[1:i-1], x[1:i-1, :], T[1:i-1], stable[1:i-1]
        end

        Function_Periodic_Orbits!(F, x[i, :], p[i], T[i], x[i-1, :], p[i-1], T[i-1], xs, ps, Ts, Δs, sol[end, :](zerosTN), dx0, n)
        Jacobian_Periodic_Orbits!(J, f!, p[i], T[i], xs, ps, Ts, sol[end,:], dx, dx0, zerosTN, n)

        j = 1

        while j < newtonite && norm(F) > newtontol
            
            try 
                Δ .= J \ F
            catch
                println("Error: El Jacobiano en la continuación es no invertible.")
                println("-> Órbitas calculadas con éxito: $(i-1)")
                return p[1:i-1, :], x[1:i-1, :], T[1:i-1], stable[1:i-1]
            end

            for k in 1:n 
                x[i,  k] -= Δ[k]
            end

            p[i] -= Δ[n+1]
            T[i] -= Δ[n+2]

            time, sol = taylorinteg(f!, x[i,:] .+ Sx, 0.0, T[i], integorder, integtol, p[i] + Sp; maxsteps = integmaxsteps)

            if length(time) > integmaxsteps
                return p[1:i-1], x[1:i-1, :], T[1:i-1], stable[1:i-1]
            end

            Function_Periodic_Orbits!(F, x[i, :], p[i], T[i], x[i-1, :], p[i-1], T[i-1], xs, ps, Ts, Δs, sol[end, :](zerosTN), dx0, n)
            Jacobian_Periodic_Orbits!(J, f!, p[i], T[i], xs, ps, Ts, sol[end,:], dx, dx0, zerosTN, n)
            
            j += 1

        end

        for i in 1:n
            for j in 1:n 
                M[i,j] = differentiate(sol[end,i], j)(zerosTN)
            end
        end
    
        if all(norm.(eigvals(M)) .<= 1.0)
            stable[i] = true
        end

        if norm(F) > maxtol
            println("La norma del sistema extendido para órbitas periódicas excedió la máxima tolerancia.")
            println("Prueba con un paso más corto")
            println("-> Órbitas calculadas con éxito: $(i-1)")
            return p[1:i-1], x[1:i-1, :], T[1:i-1], stable[1:i-1]
        end

        println("$(i): ite = $(j) : norm(F) = $(norm(F))")

        # f!(dx0, x[i,:], p[i, :], 0.0)

        # Jacobian_Periodic_Orbits!(J, f!, p[i, :], T[i], xs, ps, Ts, sol[end,:], dx, dx0, zerosTN, n)

        ΔS .= J \ V

        for j in 1:n 
            xs[j] = ΔS[j]/norm(ΔS)
        end

        ps  = ΔS[n+1]/norm(ΔS)
        Ts  = ΔS[n+2]/norm(ΔS)

        f!(dx0, x[i,:], p[i], 0.0)
        
    end

    return p, x, T, stable
    
end

function NewtonPeriodic!(i::Int64, f!::Function, F::Vector{Float64}, J::Matrix{Float64}, Φ::Matrix{Float64}, Δ::Vector{Float64}, x::Matrix{Float64}, p::Matrix{Float64}, T::Vector{Float64}, xs::Vector{Float64}, ps::Float64, Ts::Float64, Δs::Float64, dx0::Vector{Float64}, indice::Int64, n::Int64, dx::Vector{Float64}, time::Vector{Float64}, sol::Matrix{Taylor1{Float64}}, Sx::Matrix{Taylor1{Float64}}, Sp::Vector{Taylor1{Float64}}, r::Taylor1{Float64}, newtonite::Int64, newtontol::Float64, integorder::Int64, integtol::Float64, integmaxsteps::Int64)

    for j in 1:n 
        x[i,j] = x[i-1,j] + Δs*xs[j]
    end

    p[i, indice] = p[i-1, indice] + Δs*ps
    T[i] = T[i-1] + Δs*Ts

    time[2] = T[i]

    for j in 1:n
        sol .= taylorinteg(f!, x[i,:] .+ Sx[j, :], time, integorder, integtol, p[i, :] .+ r; maxsteps = integmaxsteps, parse_eqs = false)
        Φ[:, j] .= differentiate.(sol[2, :])(0.0)
    end

    sol .= taylorinteg(f!, x[i,:] .+ r, time, integorder, integtol, p[i, :] .+ Sp; maxsteps = integmaxsteps, parse_eqs = false)
    Φ[:, n+1] .= differentiate.(sol[2, :])(0.0)

    Function_Periodic_Orbits!(F, x[i, :], p[i, indice], T[i], x[i-1, :], p[i-1,indice], T[i-1], xs, ps, Ts, Δs, sol[2, :](0.0), dx0, n)
    Jacobian_Periodic_Orbits!(J, f!, p[i, :], T[i], xs, ps, Ts, Φ, sol[2, :](0.0),dx, dx0, n)

    j = 1

    while j < newtonite && norm(F) > newtontol
            
        Δ .= J \ F

        for k in 1:n 
            x[i,  k] -= Δ[k]
        end

        p[i, indice] -= Δ[n+1]
        T[i]         -= Δ[n+2]

        time[2] = T[i]

        for j in 1:n
            sol .= taylorinteg(f!, x[i,:] .+ Sx[j, :], time, integorder, integtol, p[i, :] .+ r; maxsteps = integmaxsteps, parse_eqs = false)
            Φ[:, j] .= differentiate.(sol[2, :])(0.0)
        end

        sol .= taylorinteg(f!, x[i,:] .+ r, time, integorder, integtol, p[i, :] .+ Sp; maxsteps = integmaxsteps, parse_eqs = false)
        Φ[:, n+1] .= differentiate.(sol[2, :])(0.0)

        Function_Periodic_Orbits!(F, x[i, :], p[i, indice], T[i], x[i-1, :], p[i-1,indice], T[i-1], xs, ps, Ts, Δs, sol[2, :](0.0), dx0, n)
        Jacobian_Periodic_Orbits!(J, f!, p[i, :], T[i], xs, ps, Ts, Φ, sol[2, :](0.0), dx, dx0, n)
        
        j += 1

    end

end
    

#-

function Periodic_Orbits(f!::Function, x_ini::Vector{Float64}, p_ini::Vector{Float64},t::Float64, Δs::Float64, δx::Float64, indice::Int64, N::Int64; integorder = 20, integtol = 1.e-20, newtonite = 8, newtontol = 1.e-16, maxtol = 1.e-8, taylororder = 1, integmaxsteps = 1000)

    ### Inicializando variables

    n  = length(x_ini)

    np = length(p_ini)

    max_norm = 0.0

    max_ind = 0

    dx = zeros(Float64, n)

    dx0 = zeros(Float64, n)

    stable = [false for i in 1:N]

    Dfx = zeros(Float64, n,n)

    F = zeros(Float64, n+2)

    J = zeros(Float64, n+2,n+2)

    Φ = zeros(n, n + 1)

    normF = norm(F)

    p = [p_ini[j] for i in 1:N, j in 1:np]
    x = zeros(Float64, N, n)
    T = zeros(Float64, N)

    M = zeros(n,n)

    s = Taylor1(Float64, taylororder)
    r = Taylor1(zeros(Float64, taylororder + 1))

    # S = set_variables("s", numvars = n+1, order = 1)

    # zerosTN = zeros(Float64, n+1)

    dxT1 = [s for i in 1:n]

    SxT1 = [i == j ? s : r for i in 1:n, j in 1:n]

    # R = zero(S[1])

    sol = [r for i in 1:2, j in 1:n]

    Sx = [i == j ? s : r for i in 1:n, j in 1:n]#S[1:n]
    Sp = [i == indice ? s : r for i in 1:np]

    Δ = zeros(Float64, n + 2)

    xs = zeros(Float64, n)
    ps = 0.0
    Ts = 0.0

    ΔS = zeros(Float64, n + 2)

    normΔS = 0.0

    V = [zeros(Float64, n+1) ; 1.0]

    time = zeros(2)

    ### Cálculo de la frecuencia y periodo inicial

    JacobianT1!(Dfx, f!, x_ini, p_ini, t, dxT1, SxT1, r, n)

    ω_ini = abs(imag(eigvals(Dfx)[end]))

    T_ini = 2*π/ω_ini

    time[2] = T_ini

    ### Cálculo de la órbita en el punto de equilibrio

    W = nullspace([-ω_ini*I(n) Dfx; Dfx ω_ini*I(n)])[:,1]

    ws = W[1:n]
    wc = W[n+1:end]

    # ϕ(t)  =        ws*sin(2*π*t) + wc*cos(2*π*t)
    # dϕ(t) = ω_ini*(ws*cos(2*π*t) - wc*sin(2*π*t))

    ### Asignando puntos iniciales

    x[1,:] .= x_ini
    T[1]    = T_ini

    xs .= wc / norm(wc)

    dx0 .=  δx*ω_ini*ws

    ### Calculando la primer órbita

    i = 2
    
    try 
        NewtonPeriodic!(i, f!, F, J, Φ, Δ, x, p, T, xs, ps, Ts, δx, dx0, indice, n, dx, time, sol, Sx, Sp, r, newtonite, newtontol, integorder, integtol, integmaxsteps)
    catch
        println("Error: El jacobiano del sistema se volvió singular.")
        Periodic_Data(i-1, max_norm, max_ind)
        return p[1:i-1, :], x[1:i-1, :], T[1:i-1], stable[1:i-1]
    end

    println("$(i): norm(F) = $(norm(F))")

    if norm(F) > maxtol
        println("La norma del sistema extendido para órbitas periódicas excedió la máxima tolerancia.")
        println("┌ Datos de las ramas periódicas")
        println("└ Número de puntos: $(N)")
        return p[1:i-1, :], x[1:i-1, :], T[1:i-1], stable[1]
    end

    if norm(F) > max_norm
        max_norm = normF
        max_ind = i
    end

    ΔS .= J \ V

    normΔS = norm(ΔS)

    for j in 1:n 
        xs[j] = ΔS[j]/normΔS
    end

    ps  = ΔS[n+1]/normΔS
    Ts  = ΔS[n+2]/normΔS

    f!(dx0, x[i,:], p[i, :], 0.0)

    ### Calculando el resto de puntos

    for i in 3:N
            Δs_adap = Δs
            Function_Periodic_Orbits!(F, x[i-1, :], p[i-1, indice], T[i], x[i-1, :], p[i-1,indice], T[i-1], xs, ps, Ts, Δs_adap, sol[2, :](0.0), dx0, n)
            k = 1
            while norm(F) > newtontol && k <= 10
                try 
                    NewtonPeriodic!(i, f!, F, J, Φ, Δ, x, p, T, xs, ps, Ts, Δs_adap, dx0, indice, n, dx, time, sol, Sx, Sp, r, newtonite, newtontol, integorder, integtol, integmaxsteps)
                catch
                    # println("Error: El jacobiano del sistema se volvió singular.")
                    # Periodic_Data(i-1, max_norm, max_ind)
                    # return p[1:i-1, :], x[1:i-1, :], T[1:i-1], stable[1:i-1]
                    Function_Periodic_Orbits!(F, x[i-1, :], p[i-1, indice], T[i], x[i-1, :], p[i-1,indice], T[i-1], xs, ps, Ts, Δs, sol[2, :](0.0), dx0, n)
                end
                println("$(i): norm(F) = $(norm(F))")
                Δs_adap = Δs_adap*1.e-1
                k +=1
            end
        

        if norm(F) > maxtol
            println("La norma del sistema extendido para órbitas periódicas excedió la máxima tolerancia.")
            Periodic_Data(i-1, max_norm, max_ind)
            return p[1:i-1, :], x[1:i-1, :], T[1:i-1], stable[1:i-1]
        end

        if norm(F) > max_norm
            max_norm = normF
            max_ind = i
        end

        ΔS .= J \ V
        normΔS = norm(ΔS)
        for j in 1:n 
            xs[j] = ΔS[j]/normΔS
        end
        ps  = ΔS[n+1]/normΔS
        Ts  = ΔS[n+2]/normΔS

        f!(dx0, x[i,:], p[i, :], 0.0)
        
    end

    Periodic_Data(N, max_norm, max_ind)

    return p, x, T, stable
    
end

#-