#module ImplicitFunction

    #export Implicit_function
include("nonlinear_system.jl")

#using TaylorSeries, .NonlinearSystem

function Implicit_function(f::Function, x::Float64, p::Float64, orden::Int64)
    return -derivative(f(x,p+Taylor1(orden)))(0.0)/derivative(f(x+Taylor1(orden),p))(0.0)
end

function Implicit_function(f::Function, x::Float64, p::Vector{Float64}, orden::Int64,indice::Int64)
    T = [i == indice ? Taylor1(orden) : Taylor1(0) for i in 1:length(p)]
    return -derivative(f(x,p+T))(0.0)/derivative(f(x+Taylor1(orden),p))(0.0)
end

function Implicit_function(f!::Function, x::Vector{Float64}, p::Float64, orden::Int64)
    return -inv(Jacobian(f!,x,p,orden))*Jacobian(f!,x,p+Taylor1(orden))
end

function Implicit_function(f!::Function, x::Vector{Float64}, p::Vector{Float64}, orden::Int64,indice::Int64)
    return -inv(Jacobian(f!,x,p,orden))*Jacobian(f!,x,p,orden,indice)
end

#end
end
