


function Initial_Periodic_Orbits(x_ini, x_end, x0, T0, Δs, ws, wc)
    
    pc = (x_ini - (x0 + wc*Δs)) ⋅ (2*π*ws/T0)

    psa = (x_ini - x0) ⋅ wc - Δs

    return [x_end - x_ini; pc; psa]
    
end 

#-

function Initial_Jacobian_Orbits(f!, x_ini, x_end, p, T, T0, ws, wc)
    
    n = length(x_ini)

    J = zeros(n+2,n+2)

    for i in 1:n
        J[1:n,i] .= differentiate.(x_end - x_ini, i)(zeros(n+1))
    end

    J[end-1,1:n] = transpose((2*π*ws/T0))

    J[end,1:n] = transpose(wc)

    J[1:n,end-1] .= differentiate.(x_end - x_ini, n+1)(zeros(n+1))
    
    dx = zeros(n)

    f!(dx, x_end(zeros(n+1)), p, T)

    J[1:n,end] = transpose(dx)

    return J
    
end

#-

function Initial_Jacobian_Orbits(f!, x_ini, x_end, p, T, T0, ws, wc, p_ini, indice)
    
    n = length(x_ini)

    J = zeros(n+2,n+2)

    for i in 1:n
        J[1:n,i] .= differentiate.(x_end - x_ini, i)(zeros(n+1))
    end

    J[end-1,1:n] = transpose((2*π*ws/T0))

    J[end,1:n] = transpose(wc)

    J[1:n,end-1] .= differentiate.(x_end - x_ini, n+1)(zeros(n+1))
    
    dx = zeros(n)

    f!(dx, x_end(zeros(n+1)), [i != indice ? p_ini[i] : p for i in 1:length(p_ini)], T)

    J[1:n,end] = transpose(dx)

    return J
    
end

#-

function Function_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, T0, xs, ps, Ts, Δs)
    
    dx_old = zeros(length(x0))
    
    f!(dx_old,x0, p, 0.0)
    
    pc = (x_ini - x0) ⋅ dx_old

    psa = (x_ini - x0) ⋅ xs + (p - p0)*ps + (T - T0)*Ts - Δs
    
    return [x_end - x_ini ; pc; psa]
    
end

#-

function Function_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, T0, xs, ps, Ts, Δs, p_ini, indice)
    
    dx_old = zeros(length(x0))
    
    f!(dx_old,x0,[i != indice ? p_ini[i] : p0 for i in 1:length(p_ini)],0.0)
    
    pc = (x_ini - x0) ⋅ dx_old

    psa = (x_ini - x0) ⋅ xs + (p - p0)*ps + (T - T0)*Ts - Δs
    
    return [x_end - x_ini ; pc; psa]
    
end

#-

function Jacobian_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, xs, ps, Ts)
    
    n = length(x0)

    dx_old = zeros(n)
    
    f!(dx_old,x0, p0, 0.0)

    J = zeros(n+2,n+2)

    for i in 1:n
        J[1:n,i] = differentiate.(x_end - x_ini, i)(zeros(n+1))
    end

    J[end - 1,1:n] = transpose(dx_old)

    J[end,1:n] = transpose(xs)

    J[1:n,end-1] = differentiate.(x_end - x_ini, n+1)(zeros(n+1))

    J[end,end-1] = ps
    
    dx = zeros(n)

    f!(dx, x_end(zeros(n+1)), p, T)

    J[1:n,end] = transpose(dx)

    J[end,end] = Ts

    return J
    
end

#-

function Jacobian_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, xs, ps, Ts, p_ini, indice)
    
    n = length(x0)

    dx_old = zeros(n)
    
    f!(dx_old,x0,[i != indice ? p_ini[i] : p0 for i in 1:length(p_ini)],0.0)

    J = zeros(n+2,n+2)

    for i in 1:n
        J[1:n,i] = differentiate.(x_end - x_ini, i)(zeros(n+1))
    end

    J[end - 1,1:n] = transpose(dx_old)

    J[end,1:n] = transpose(xs)

    J[1:n,end-1] = differentiate.(x_end - x_ini, n+1)(zeros(n+1))

    J[end,end-1] = ps
    
    dx = zeros(n)

    f!(dx, x_end(zeros(n+1)), [i != indice ? p_ini[i] : p for i in 1:length(p_ini)], T)

    J[1:n,end] = transpose(dx)

    J[end,end] = Ts

    return J
    
end

#-

function Periodic_Orbits(f!, x_bif, p_bif, t, Δs, p_ini, p_fin, indice; N = 100, integorder = 20, integtol = 1.e-20, newtonite = 8, newtontol = 1.e-16)

    n = length(x_bif)

    Jx = JacobianT1(f!,x_bif, p_bif, t, p_ini, indice)

    ω = abs(imag(eigvals(Jx)[end]))

    T_bif = 2*π/ω

    X  = [x_bif]
    P  = [p_bif]
    TP = [T_bif]

    W = nullspace([-ω*I(n) Jx; Jx ω*I(n)])

    ws = W[:,1][1:n]
    wc = W[:,1][n+1:end]

    x0 = x_bif
    p0 = p_bif
    T0 = T_bif

    x = x_bif + Δs*wc
    p = p_bif
    T = T_bif

    i = 1

    S = set_variables("s", numvars = n + 1, order = 1)

    time, sol = taylorinteg(f!, x .+ S[1:n], 0.0, T, integorder, integtol, [k != indice ? p_ini[k] : p + S[end] for k in 1:length(p_ini)])

    x_ini = sol[1,:]
    x_end = sol[end,:]

    IPO = Initial_Periodic_Orbits(x_ini(zeros(n+1)), x_end(zeros(n+1)), x0, T0, Δs, ws, wc)
    IJO = Initial_Jacobian_Orbits(f!, x_ini, x_end, p, T, T0, ws, wc, p_ini, indice)

    while i <= newtonite && norm(IPO) > newtontol

        Δ = -inv(IJO)*IPO

        x += Δ[1:n]
        p += Δ[end-1]
        T += Δ[end]
        
        time, sol = taylorinteg(f!, x .+ S[1:n], 0.0, T, integorder, integtol, [k != indice ? p_ini[k] : p + S[end] for k in 1:length(p_ini)])

        x_ini = sol[1,:]
        x_end = sol[end,:]

        IPO = Initial_Periodic_Orbits(x_ini(zeros(n+1)), x_end(zeros(n+1)), x0, T0, Δs, ws, wc)
        IJO = Initial_Jacobian_Orbits(f!, x_ini, x_end, p, T, T0, ws, wc, p_ini, indice)

        i += 1

    end

    x0 = deepcopy(x)
    p0 = deepcopy(p)
    T0 = deepcopy(T)

    xs = zeros(n)
    ps = 0.0
    Ts = 0.0

    ΔS = nullspace(Jacobian_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, xs, ps, Ts, p_ini, indice)[1:end-1,:])

    xs = ΔS[1:n]
    ps = ΔS[end-1]
    Ts = ΔS[end]

    j = 1

    while p_ini[indice] <= p <= p_fin && j <= N

        FPO = Function_Periodic_Orbits(f!, x_ini(zeros(n+1)), x_end(zeros(n+1)), p, T, x0, p0, T0, xs, ps, Ts, Δs, p_ini, indice)
        JPO = Jacobian_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, xs, ps, Ts, p_ini, indice)

        i = 1
        
        while i <= newtonite && norm(FPO) > newtontol
            
            Δ = -inv(JPO)*FPO

            x += Δ[1:n]
            p += Δ[end-1]
            T += Δ[end]

            time, sol = taylorinteg(f!, x .+ S[1:n], 0.0, T, integorder, integtol, [k != indice ? p_ini[k] : p + S[end] for k in 1:length(p_ini)])

            x_ini = sol[1,:]
            x_end = sol[end,:]
            
            FPO = Function_Periodic_Orbits(f!, x_ini(zeros(n+1)), x_end(zeros(n+1)), p, T, x0, p0, T0, xs, ps, Ts, Δs, p_ini, indice)
            JPO = Jacobian_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, xs, ps, Ts, p_ini, indice)
            
            i += 1
        end

        x0 = deepcopy(x)
        p0 = deepcopy(p)
        T0 = deepcopy(T)
    
        ΔS = nullspace(JPO[1:end-1,:])
    
        xs = ΔS[1:n]
        ps = ΔS[end-1]
        Ts = ΔS[end]

        push!(X,x)
        push!(P,p)
        push!(TP,T)

        println(j)
        
        j += 1
        
    end

    return P, X, TP
    
end

#-

function Periodic_Orbits(f!, x_bif, p_bif, t, Δs, p_ini, p_fin, indice; N = 100, integorder = 20, integtol = 1.e-20, newtonite = 8, newtontol = 1.e-16)

    n = length(x_bif)

    Jx = JacobianT1(f!,x_bif, p_bif, t, p_ini, indice)

    ω = abs(imag(eigvals(Jx)[end]))

    T_bif = 2*π/ω

    X  = [x_bif]
    P  = [p_bif]
    TP = [T_bif]

    W = nullspace([-ω*I(n) Jx; Jx ω*I(n)])

    ws = W[:,1][1:n]
    wc = W[:,1][n+1:end]

    x0 = x_bif
    p0 = p_bif
    T0 = T_bif

    x = x_bif + Δs*wc
    p = p_bif
    T = T_bif

    i = 1

    S = set_variables("s", numvars = n + 1, order = 1)

    time, sol = taylorinteg(f!, x .+ S[1:n], 0.0, T, integorder, integtol, [k != indice ? p_ini[k] : p + S[end] for k in 1:length(p_ini)])

    x_ini = sol[1,:]
    x_end = sol[end,:]

    IPO = Initial_Periodic_Orbits(x_ini(zeros(n+1)), x_end(zeros(n+1)), x0, T0, Δs, ws, wc)
    IJO = Initial_Jacobian_Orbits(f!, x_ini, x_end, p, T, T0, ws, wc, p_ini, indice)

    while i <= newtonite && norm(IPO) > newtontol

        Δ = -inv(IJO)*IPO

        x += Δ[1:n]
        p += Δ[end-1]
        T += Δ[end]
        
        time, sol = taylorinteg(f!, x .+ S[1:n], 0.0, T, integorder, integtol, [k != indice ? p_ini[k] : p + S[end] for k in 1:length(p_ini)])

        x_ini = sol[1,:]
        x_end = sol[end,:]

        IPO = Initial_Periodic_Orbits(x_ini(zeros(n+1)), x_end(zeros(n+1)), x0, T0, Δs, ws, wc)
        IJO = Initial_Jacobian_Orbits(f!, x_ini, x_end, p, T, T0, ws, wc, p_ini, indice)

        i += 1

    end

    x0 = deepcopy(x)
    p0 = deepcopy(p)
    T0 = deepcopy(T)

    xs = zeros(n)
    ps = 0.0
    Ts = 0.0

    ΔS = nullspace(Jacobian_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, xs, ps, Ts, p_ini, indice)[1:end-1,:])

    xs = ΔS[1:n]
    ps = ΔS[end-1]
    Ts = ΔS[end]

    j = 1

    while p_ini[indice] <= p <= p_fin && j <= N

        FPO = Function_Periodic_Orbits(f!, x_ini(zeros(n+1)), x_end(zeros(n+1)), p, T, x0, p0, T0, xs, ps, Ts, Δs, p_ini, indice)
        JPO = Jacobian_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, xs, ps, Ts, p_ini, indice)

        i = 1
        
        while i <= newtonite && norm(FPO) > newtontol
            
            Δ = -inv(JPO)*FPO

            x += Δ[1:n]
            p += Δ[end-1]
            T += Δ[end]

            time, sol = taylorinteg(f!, x .+ S[1:n], 0.0, T, integorder, integtol, [k != indice ? p_ini[k] : p + S[end] for k in 1:length(p_ini)])

            x_ini = sol[1,:]
            x_end = sol[end,:]
            
            FPO = Function_Periodic_Orbits(f!, x_ini(zeros(n+1)), x_end(zeros(n+1)), p, T, x0, p0, T0, xs, ps, Ts, Δs, p_ini, indice)
            JPO = Jacobian_Periodic_Orbits(f!, x_ini, x_end, p, T, x0, p0, xs, ps, Ts, p_ini, indice)
            
            i += 1
        end

        x0 = deepcopy(x)
        p0 = deepcopy(p)
        T0 = deepcopy(T)
    
        ΔS = nullspace(JPO[1:end-1,:])
    
        xs = ΔS[1:n]
        ps = ΔS[end-1]
        Ts = ΔS[end]

        push!(X,x)
        push!(P,p)
        push!(TP,T)

        println(j)
        
        j += 1
        
    end

    return P, X, TP
    
end