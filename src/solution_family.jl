module SolutionFamily

    export Solution_family

    include("psarc.jl")
    include("newton.jl")

    import .PArcLength, .TaylorNewton
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
end