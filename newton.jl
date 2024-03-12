module TaylorNewton

export Newton

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

function Newton(f::Function,x0::Float64,p0::Vector{Taylor1})
    p_new = p0
    i = 1
    while i <= 30 && abs(f(x0,p_new(0.0))) > 1.e-16
        p_old = p_new
        p_new = p_old - (f(x0,p_old)/derivative(f(x0,p_old)))(0.0)
        i+=1
    end
    return p_new(0.0)
end

end