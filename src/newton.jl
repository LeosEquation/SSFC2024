module TaylorNewton

    export Newton

    include("nonlinear_system.jl")

    using TaylorSeries, LinearAlgebra, .NonlinearSystem

    function Newton(f::Function,x0::Taylor1,p0::Union{Float64,Vector{Float64}})
        x_new = x0
        i = 1
        while i <= 30 && abs(f(x_new(0.0),p0)) > 1.e-16
            x_old = x_new
            x_new = x_old - (f(x_old,p0)/derivative(f(x_old,p0)))(0.0)
            i+=1
        end
        return x_new(0.0)
    end

    function Newton(f::Function,x0::Float64,p0::Taylor1)
        p_new = p0
        i = 1
        while i <= 30 && abs(f(x0,p_new(0.0))) > 1.e-16
            p_old = p_new
            p_new = p_old - (f(x0,p_old)/derivative(f(x0,p_old)))(0.0)
            i+=1
        end
        return p_new(0.0)
    end

    function Newton(f::Function,x0::Float64,p0::Vector{Taylor1{Float64}},indice::Int64)
        p_new = p0
        i = 1
        while i <= 30 && abs(f(x0,p_new(0.0))) > 1.e-16
            p_old = p_new
            p_new = p_old - [i == indice ? (f(x0,p_old)/derivative(f(x0,p_old)))(0.0) : Taylor1(0) for i in 1:length(p0)]
            i+=1
        end
        return p_new(0.0)
    end

    function Newton(f!::Function,x0::Vector{Float64},p0::Union{Float64,Vector{Float64}},orden::Int64)
        x_new = x0
        i = 1
        dx_new = [Taylor(0) for i in 1:length(x0)]
        f!(dx,x0,p0)
        while i <= 30 && norm(dx_new(0.0)) > 1.e-16
            x_old = x_new
            x_new = x_old - inv(NonlinearSystem.Jacobian(f!,x_old(0.0),p0,orden))*dx_new(0.0)
            i+=1
            f!(dx_new,x_new,p0)
        end
        return x_new(0.0)
    end

    function Newton(f!::Function,x0::Vector{Float64},p0::Taylor1)
        p_new = p0
        i = 1
        dx_new = [Taylor(0) for i in 1:length(x0)]
        f!(dx,x0,p0)
        while i <= 30 && norm(dx_new(0.0)) > 1.e-16
            p_old = p_new
            p_new = p_old - inv(NonlinearSystem.Jacobian(f!,x0,p_old))*dx_new(0.0)
            i+=1
        end
        return p_new(0.0)
    end

    function Newton(f!::Function,x0::Float64,p0::Vector{Taylor1{Float64}},orden::Int64,indice::Int64)
        p_new = p0
        i = 1
        dx_new = [Taylor(0) for i in 1:length(x0)]
        f!(dx,x0,p0)
        while i <= 30 && norm(dx_new(0.0)) > 1.e-16
            p_old = p_new
            p_new = p_old - [j == indice ? inv(NonlinearSystem.Jacobian(f!,x0,p_old,orden,indice))*dx_new(0.0) : 0.0 for j in 1:length(p)]
            i+=1
        end
        return p_new(0.0)
    end


end