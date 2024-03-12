module ImplicitFunction

export Implicit_function

using TaylorSeries

function Implicit_function(f::Function, x::Float64, p::Float64, orden::Int64, indice::Int64)
    return -derivative(f(x,p+Taylor1(orden)))(0.0)/derivative(f(x+Taylor1(orden),p))(0.0)
end

end