module PArcLength

    export first_step, step

    include("implicit_function.jl")
    include("nonlinear_system.jl")

    import .ImplicitFunction, .NonlinearSystem

    using TaylorSeries, LinearAlgebra


    function first_step(f::Function,x_ini::Float64,p_ini::Float64,p_fin::Float64,orden::Int64)

        x_s = 1/sqrt(1/(ImplicitFunction.Implicit_function(f,x_ini,p_ini,orden)^2) + 1)
        p_s = sign(p_fin - p_ini)/sqrt(ImplicitFunction.Implicit_function(f,x_ini,p_ini,orden)^2 + 1)

        return x_s, p_s

    end

    function step(f::Function,x_ini::Float64, p_ini::Float64, x_s::Float64,p_s::Float64,orden::Int64)
        t = Taylor1(orden)
        
        f_x = derivative(f(x_ini+t,p_ini))(0.0)
        f_p = derivative(f(x_ini,p_ini+t))(0.0)

        A = [f_x f_p;
            x_s p_s]

        b = transpose([0 1])

        x_s_new, p_s_new = A\b

        return x_s_new, p_s_new
    end



    function first_step(f::Function,x_ini::Float64,p_ini::Vector{Float64},p_fin::Float64,orden::Int64, indice::Int64)
        x_s = 1/sqrt(1/(ImplicitFunction.Implicit_function(f,x_ini,p_ini,orden,indice)^2) + 1)
        p_s = sign(p_fin - p_ini[indice])/sqrt(ImplicitFunction.Implicit_function(f,x_ini,p_ini,orden,indice)^2 + 1)

        return x_s, p_s

    end

    function step(f::Function,x_ini::Float64, p_ini::Vector{Float64}, x_s::Float64,p_s::Float64,orden::Int64,indice::Int64)
        T = [i == indice ? Taylor1(orden) : Taylor1(0) for i in 1:length(p_ini)]

        f_x = derivative(f(x_ini+Taylor1(orden),p_ini))(0.0)
        f_p = derivative(f(x_ini,p_ini+T))(0.0)

        A = [f_x f_p;
            x_s p_s]

        b = transpose([0 1])

        x_s_new, p_s_new = A\b

        return x_s_new, p_s_new
    end




    function first_step(f!::Function,x_ini::Vector{Float64},p_ini::Float64,p_fin::Float64,orden::Int64)

        p_s = sign(p_fin - p_ini)/sqrt(norm(ImplicitFunction.Implicit_function(f!,x_ini,p_ini,orden))^2 + 1)
        x_s = abs(p_s)*ImplicitFunction.Implicit_function(f!,x_ini,p_ini,orden)
        return x_s, p_s

    end

    function step(f!::Function,x_ini::Vector{Float64}, p_ini::Float64, x_s::Vector{Float64},p_s::Float64,orden::Int64)
        t = Taylor1(orden)
        
        f_x = NonlinearSystem.Jacobian(f!,x_ini,p_ini,orden)
        f_p = NonlinearSystem.Jacobian(f!,x_ini,p_ini+Taylor1(orden))

        A = [f_x f_p;
            transpose(x_s) p_s]

        m = length(x_s)+1

        b = [i != m ? 0.0 : 1.0 for i in 1:m]

        x_s_new, p_s_new = A\b

        return x_s_new, p_s_new
    end



    function first_step(f!::Function,x_ini::Vector{Float64},p_ini::Vector{Float64},p_fin::Float64,orden::Int64, indice::Int64)
        p_s = sign(p_fin - p_ini[indice])/sqrt(ImplicitFunction.Implicit_function(f!,x_ini,p_ini,orden,indice)^2 + 1)
        x_s = abs(p_s)*ImplicitFunction.Implicit_function(f!,x_ini,p_ini,orde,indice)
        return x_s, p_s

    end

    function step(f!::Function,x_ini::Vector{Float64}, p_ini::Vector{Float64}, x_s::Vector{Float64},p_s::Float64,orden::Int64,indice::Int64)
        T = [i == indice ? Taylor1(orden) : Taylor1(0) for i in 1:length(p_ini)]

        f_x = NonlinearSystem.Jacobian(f!,x_ini,p_ini,orden)
        f_p = NonlinearSystem.Jacobian(f!,x_ini,p_ini,orden,indice)

        A = [f_x f_p;
            transpose(x_s) p_s]

        m = length(x_s)+1

        b = [i != m ? 0.0 : 1.0 for i in 1:m]

        x_s_new, p_s_new = A\b

        return x_s_new, p_s_new
    end


end