


function Periodic_Data(N, max_norm, max_ind)
    println("┌ Datos de las ramas periódicas")

    println("├ Número de puntos: $(N)")

    println("└ Máximo valor de la integral ∫₀ᵀF(x,p)dt : $(max_norm) (punto # $(max_ind))")
end


function Periodic_Orbits(f!::Function, x_ini::Vector{Float64}, p_ini::Float64, t::Float64, Δs::Float64, δx::Float64, N::Int64; integorder = 20, integtol = 1.e-20, newtonite = 8, newtontol = 1.e-16, maxtol = 1.e-8, taylororder = 1, integmaxsteps = 1000)

    ### Inicializando variables

    n  = length(x_ini)

    max_norm = 0.0

    max_ind = 0

    dx = zeros(Float64, n)

    dx0 = zeros(Float64, n)

    stable = [false for i in 1:N]

    Dfx = zeros(Float64, n,n)

    F = zeros(Float64, n+2)

    J = zeros(Float64, n+2,n+2)

    normF = norm(F)

    p = zeros(Float64, N)
    x = zeros(Float64, N, n)
    T = zeros(Float64, N)

    M = zeros(n,n)

    s = Taylor1(Float64, taylororder)
    r = Taylor1(zeros(Float64, taylororder + 1))

    S = set_variables("s", numvars = n+1, order = 1)

    zerosTN = zeros(Float64, n+1)

    dxT1 = [s for i in 1:n]

    SxT1 = [i == j ? s : r for i in 1:n, j in 1:n]

    R = zero(S[1])

    sol = [R for i in 1:2, j in 1:n]

    Sx = S[1:n]
    Sp = S[n+1]

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
    p[1]    = p_ini

    xs .= wc / norm(wc)

    dx0 .=  δx*ω_ini*ws

    ### Calculando la primer órbita

    i = 2
    
    try 
        NewtonPeriodic!(i, f!, F, J, Δ, x, p, T, xs, ps, Ts, δx, dx0, n, dx, time, sol, Sx, Sp, zerosTN, newtonite, newtontol, integorder, integtol, integmaxsteps)
    catch
        println("Error: El jacobiano del sistema se volvió singular.")
        Periodic_Data(i-1, max_norm, max_ind)
        return p[1:i-1], x[1:i-1, :], T[1:i-1], stable[1:i-1]
    end

    println("$(i): norm(F) = $(norm(F))")

    if norm(F) > maxtol
        println("La norma del sistema extendido para órbitas periódicas excedió la máxima tolerancia.")
        println("┌ Datos de las ramas periódicas")
        println("└ Número de puntos: $(N)")
        return p[1:i-1], x[1:i-1, :], T[1:i-1], stable[1]
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

    f!(dx0, x[i,:], p[i], 0.0)

    ### Calculando el resto de puntos

    for i in 3:N

        try 
            Δs_adap = Δs
            Function_Periodic_Orbits!(F, x[i, :], p[i], T[i], x[i-1, :], p[i-1], T[i-1], xs, ps, Ts, Δs_adap, sol[2, :](zerosTN), dx0, n)
            k = 1
            while (norm(F) > newtontol || isnan(norm(F)))  && k <= 10
                NewtonPeriodic!(i, f!, F, J, Δ, x, p, T, xs, ps, Ts, Δs_adap, dx0, n, dx, time, sol, Sx, Sp, zerosTN, newtonite, newtontol, integorder, integtol, integmaxsteps)
                println("$(i): norm(F) = $(norm(F))")
                Δs_adap = Δs_adap*1.e-1
                k +=1
            end
        catch
            println("Error: El jacobiano del sistema se volvió singular.")
            Periodic_Data(i-1, max_norm, max_ind)
            return p[1:i-1], x[1:i-1, :], T[1:i-1], stable[1:i-1]
        end

        if norm(F) > maxtol || isnan(norm(F)) || isinf(norm(F))
            println("La norma del sistema extendido para órbitas periódicas excedió la máxima tolerancia.")
            Periodic_Data(i-1, max_norm, max_ind)
            return p[1:i-1], x[1:i-1, :], T[1:i-1], stable[1:i-1]
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
        f!(dx0, x[i,:], p[i], 0.0)
        
    end

    Periodic_Data(N, max_norm, max_ind)

    return p, x, T, stable
    
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

    normF = norm(F)

    p = [p_ini[j] for i in 1:N, j in 1:np]
    x = zeros(Float64, N, n)
    T = zeros(Float64, N)

    M = zeros(n,n)

    s = Taylor1(Float64, taylororder)
    r = Taylor1(zeros(Float64, taylororder + 1))

    S = set_variables("s", numvars = n+1, order = 1)

    zerosTN = zeros(Float64, n+1)

    dxT1 = [s for i in 1:n]

    SxT1 = [i == j ? s : r for i in 1:n, j in 1:n]

    R = zero(S[1])

    sol = [R for i in 1:2, j in 1:n]

    Sx = S[1:n]
    Sp = [i == indice ? S[n+1] : R for i in 1:np]

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
        NewtonPeriodic!(i, f!, F, J, Δ, x, p, T, xs, ps, Ts, δx, dx0, indice, n, dx, time, sol, Sx, Sp, zerosTN, newtonite, newtontol, integorder, integtol, integmaxsteps)
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

        try 
            Δs_adap = Δs
            Function_Periodic_Orbits!(F, x[i, :], p[i, indice], T[i], x[i-1, :], p[i-1,indice], T[i-1], xs, ps, Ts, Δs_adap, sol[2, :](zerosTN), dx0, n)
            k = 1
            while (norm(F) > newtontol || isnan(norm(F)))  && k <= 10
                NewtonPeriodic!(i, f!, F, J, Δ, x, p, T, xs, ps, Ts, Δs_adap, dx0, indice, n, dx, time, sol, Sx, Sp, zerosTN, newtonite, newtontol, integorder, integtol, integmaxsteps)
                println("$(i): norm(F) = $(norm(F))")
                Δs_adap = Δs_adap*1.e-1
                k +=1
            end
        catch
            println("Error: El jacobiano del sistema se volvió singular.")
            Periodic_Data(i-1, max_norm, max_ind)
            return p[1:i-1, :], x[1:i-1, :], T[1:i-1], stable[1:i-1]
        end

        if norm(F) > maxtol || isnan(norm(F)) || isinf(norm(F))
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