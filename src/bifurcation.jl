include("newton.jl")
include("implicit_function.jl")

using TaylorSeries

function Bifurcation_point(f::Function,x0::Float64,p0::Float64,orden::Int64)
    p_new = p0
    x_new = x0
    x, p = set_variables("x p",order = orden)
    i = 1
    while i <= 30 && abs(Implicit_function(f,x_new,p_new,orden)) > 1.e-16
        p_old = p_new
        x_old = x_new
        x_new = Newton(f,x_old + Taylor1(orden),p_old)
        p_new = p_old - derivative(f(x_new,p_old+Taylor1(orden)))(0.0)*differentiate(f(x_new+x,p_old+p), (1,1))(0.0,0.0) / (derivative(f(x_new + Taylor1(orden),p_old))(0.0)*differentiate(f(x_new,p_old+Taylor1(orden)), 2))(0.0)
        i+=1
    end
    return x_new
end