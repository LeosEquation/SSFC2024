module PArcLength

export SolFam

using TaylorSeries, LinearAlgebra

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

function Newton(f::Function,x0::Vector{TaylorN},p0::Union{Float64,Vector{Float64}})
    x_new = x0
    N = length(x0)
    i = 1
    while i <= 30 && abs(f(x_new(zeros(N)),p0)) > 1.e-16
        x_old = x_new
        x_new = x_old - inv(jacobian(x_old)) * f(x_old(zeros(N)),p0)
        i+=1
    end
    return x_new(zeros(N))
end

function Newton(f::Function,x0::Union{Float64,Vector{Float64}},p0::Taylor1)
    p_new = p0
    i = 1
    while i <= 30 && abs(f(x0,p_new(0.0))) > 1.e-16
        p_old = p_new
        p_new = p_old - (f(x0,p_old)/derivative(f(x0,p_old)))(0.0)
        i+=1
    end
    return p_new(0.0)
end

x_p(f::Function, x::Float64, p::Float64, orden::Int64) = -derivative(f(x,p+Taylor1(orden)))(0.0)/derivative(f(x+Taylor1(orden),p))(0.0)


"""function IFT(f::Function, x::Float64, p::Union{Float64,Vector{Float64}}, orden::Int64, indice::Int64)
    
    if length(p) > 1
        T = [Taylor1(orden):Taylor1(1)(0.0) ? i = indice for i in 1:length(p)]
    else 
        T = Taylor1(orden)
    end

    return -derivative(f(x,p+T))(0.0)/derivative(f(x+Taylor1(orden),p))(0.0)
end"""



function first_step(f::Function,x_ini::Float64,p_ini::Float64,p_fin::Float64,orden::Int64)
    t = Taylor1(orden)

    x_s = 1/sqrt(1/(x_p(f,x_ini,p_ini,orden)^2) + 1)
    p_s = sign(p_fin - p_ini)/sqrt(x_p(f,x_ini,p_ini,orden)^2 + 1)

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

function SolFam(f::Function, x_ini::Float64, p_ini::Float64, Δs::Float64, p_fin::Float64,orden::Int64; tol = 1.e-16, N = 10000)
        
    X = Float64[]
    P = Float64[]

    t = Taylor1(orden)

    push!(P,p_ini)
    push!(X,x_ini)

    x_s_new, p_s_new = first_step(f, x_ini, p_ini, p_fin,orden)
    x_old = x_ini + x_s_new*Δs
    p_old = p_ini + p_s_new*Δs

    
    x_new = Newton(f,x_old + t,p_old)
    p_new = Newton(f,x_new,p_old + t)
    if abs(f(x_new,p_new)) <= 1.e-16
        push!(P,p_new)
        push!(X,x_new)
        x_old = x_new
        p_old = p_new
        x_s_old = x_s_new
        p_s_old = p_s_new
    else
        throw(ArgumentError("No se pudo calcular el primer paso, prueba con valores iniciales distintos"))
    end

    i = 0

    while p_ini <= p_old <= p_fin && i <= N
        x_s_new, p_s_new = step(f, x_old, p_old, x_s_old, p_s_old ,orden)
        x_old += x_s_new*Δs
        p_old += p_s_new*Δs
        x_new = Newton(f,x_old + t,p_old)
        p_new = Newton(f,x_new,p_old + t)
        push!(P,p_new)
        push!(X,x_new)
        x_old = x_new
        p_old = p_new
        x_s_old = x_s_new
        p_s_old = p_s_new
        i+=1  
    end

    return P[1:end-1],X[1:end-1]

end

end