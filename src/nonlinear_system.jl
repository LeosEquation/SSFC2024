#module NonlinearSystem

#    export Jacobian

#    using TaylorSeries, LinearAlgebra

function Jacobian(f!::Function, x::Vector{Float64}, p::Union{Float64,Vector{Float64}},orden::Int64)
    J = zeros(length(x),length(x))
    for i in 1:length(x)
        dx = [Taylor1(0) for i in 1:length(x)]
        f!(dx,x + [i == j ? Taylor1(orden) : Taylor1([0.0],orden) for j in 1:length(x)],p)
        for j in 1:length(x)
            J[j,i] = derivative(dx[j])(0.0)
        end
    end
    #println(J)
    return J
end

function Jacobian(f!::Function, x::Vector{Float64}, p::Taylor1)
    J = zeros(length(x))
    dx = [Taylor1(0) for i in 1:length(x)]
    f!(dx,x,p)
    for i in 1:length(x)
        J[i] = derivative(dx[i] + Taylor1([0.0],(length(p)-1)))(0.0)
    end
    #println("J_p = $(J)")
    return J
end

function Jacobian(f!::Function, x::Vector{Float64}, p::Vector{Float64},orden::Int64,indice::Int64)
    T = [i == indice ? Taylor1(orden) : Taylor1([0.0],orden) for i in 1:length(p)]
    J = zeros(length(x))
    dx = [Taylor1(0) for i in 1:length(x)]
    f!(dx,x,p + T)
    for i in 1:length(x)
        J[i] = derivative(dx[i] + Taylor1([0.0],(length(p[indice])-1)))(0.0)
    end
    return J
end

#end