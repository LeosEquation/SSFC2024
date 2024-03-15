module SolutionFamily

    export Solution_family

    include("psarc.jl")
    include("newton.jl")
    include("implicit_function.jl")

    import .PArcLength, .TaylorNewton, .ImplicitFunction
    using TaylorSeries 

    function Solution_family(f::Function, x_ini::Float64, p_ini::Float64, Δs::Float64, p_fin::Float64,orden::Int64; tol = 1.e-16, N = 10000)
            
        X = Float64[]
        P = Float64[]

        t = Taylor1(orden)

        push!(P,p_ini)
        push!(X,x_ini)

        x_s, p_s = PArcLength.first_step(f, x_ini, p_ini, p_fin,orden)
        x = x_ini + x_s*Δs
        p = p_ini + p_s*Δs

        
        x = TaylorNewton.Newton(f,x + t,p)
        p = TaylorNewton.Newton(f,x,p + t)
        if abs(f(x,p)) <= 1.e-16
            push!(P,p)
            push!(X,x)
        else
            throw(ArgumentError("No se pudo calcular el primer paso, prueba con valores iniciales distintos"))
        end
        

        i = 0

        while p_ini <= p <= p_fin && i <= N
            x_s, p_s = PArcLength.step(f, x, p, x_s, p_s ,orden)
            x = TaylorNewton.Newton(f,x + x_s*Δs + t,p + p_s*Δs)
            p = TaylorNewton.Newton(f,x,p + p_s*Δs + t)
            push!(P,p)
            push!(X,x)
            i+=1  
        end

        return P[1:end-1],X[1:end-1]

    end

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

        x_s, p_s = PArcLength.first_step(f, x_ini, p_ini, p_fin,orden,indice)
        x = x_ini + x_s*Δs
        p = p_ini + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p_ini)]

        
        x = TaylorNewton.Newton(f,x + t,p)
        p = TaylorNewton.Newton(f,x,p + T,indice)
        if abs(f(x,p)) <= 1.e-16
            push!(P,p[indice])
            push!(X,x)
        else
            throw(ArgumentError("No se pudo calcular el primer paso, prueba con valores iniciales distintos"))
        end
        

        i = 0

        while p_ini[indice] <= p[indice] <= p_fin && i <= N
            x_s, p_s = PArcLength.step(f, x, p, x_s, p_s ,orden,indice)
            x = TaylorNewton.Newton(f,x + x_s*Δs + t,p + [i == indice ? p_s*Δs : 0.0 for i in 1:length(p)])
            p = TaylorNewton.Newton(f,x,p + [i == indice ? p_s*Δs + t : Taylor1(0) for i in 1:length(p)],indice)
            push!(P,p[indice])
            push!(X,x)
            i+=1  
        end

        return P[1:end-1],X[1:end-1]

    end

end