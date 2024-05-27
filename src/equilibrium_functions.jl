

#using TaylorSeries, LinearAlgebra

function Equilibrium_initial_values(f::Function,x_ini::Float64,p_ini::Float64,t::Float64,p_fin::Float64)
    s = Taylor1(1)
    r = Taylor1([0.0,0.0],1)
    
    D = -derivative(f(x_ini+r,p_ini+s,t+r))(0.0)/derivative(f(x_ini+s,p_ini+r,t+r))(0.0)

    p_s = 1/sqrt(norm(D)^2 + 1)
    x_s = p_s*D
    return x_s, sign(p_fin - p_ini)*p_s
end

function Equilibrium_initial_values(f::Function,x_ini::Float64,p_ini::Vector{Float64},t::Float64,p_fin::Float64,indice::Int64)
    
    s = Taylor1(1)
    r = Taylor1([0.0,0.0],1)
    S = [i == indice ? s : r for i in 1:length(p)]
    D = -derivative(f(x_ini+r,p_ini+S,t+r))(0.0)/derivative(f(x_ini+s,p_ini .+ r,t+r))(0.0)

    p_s = 1/sqrt(norm(D)^2 + 1)
    x_s = p_s * D
    return x_s, sign(p_fin - p_ini[indice])*p_s
end

function Equilibrium_initial_values(f!::Function,x_ini::Vector{Float64},p_ini::Float64,t::Float64,p_fin::Float64)
    
    n = length(x_ini)
    s = Taylor1(2)
    r = Taylor1([0.0,0.0],2)

    Jx = zeros(n,n)
    for i in 1:n
        for j in 1:n
            dx = [s for i in 1:n]
            f!(dx,x_ini + [k == j ? s : r for k in 1:n], p_ini + r, t + r)
            Jx[i,j] = differentiate(dx[i])(0.0)
        end
    end

    Jp = zeros(n)
    for i in 1:n
        dx = [s for i in 1:n]
        f!(dx,x_ini .+ r, p_ini + s, t + r)
        Jp[i] = differentiate(dx[i])(0.0)
    end

    D = - inv(Jx) * Jp

    p_s = 1/sqrt(norm(D)^2 + 1)
    x_s = p_s*D

    return x_s, sign(p_fin - p_ini)*p_s
end

function Equilibrium_initial_values(f!::Function,x_ini::Vector{Float64},p_ini::Vector{Float64},t::Float64,p_fin::Float64,indice::Int64)
    
    n = length(x_ini)
    s = Taylor1(2)
    r = Taylor1([0.0,0.0],2)
    
    Jx = zeros(n,n)
    for i in 1:n
        for j in 1:n
            dx = [s for i in 1:n]
            f!(dx,x_ini + [k == j ? s : r for k in 1:n], p_ini .+ r, t + r)
            Jx[i,j] = differentiate(dx[i])(0.0)
        end
    end

    Jp = zeros(n)
    for i in 1:n
        dx = [s for i in 1:n]
        f!(dx,x_ini .+ r, p_ini + [k == indice ? s : r for k in 1:length(p_ini)], t + r)
        Jp[i] = differentiate(dx[i])(0.0)
    end

    D = - inv(Jx) * Jp

    p_s = 1/sqrt(norm(D)^2 + 1)
    x_s = p_s*D

    return x_s, sign(p_fin - p_ini[indice])*p_s
end


function Equilibrium_Function(f::Function,x::Float64,p::Float64,t::Float64,x0::Float64,p0::Float64,x0_s::Float64,p0_s::Float64,Δs::Float64)
    
    EF1 = f(x,p,t)
    EF2 = (x-x0)*x0_s + (p-p0)*p0_s - Δs

    return [EF1 ; EF2]
end

function Equilibrium_Function(f::Function,x::Float64,p::Float64,t::Float64,x0::Float64,p0::Float64,x0_s::Float64,p0_s::Float64,Δs::Float64,p_ini::Vector{Float64},indice::Int64)
    

    EF1 = f(x,[i == indice ? p0 : p_ini[i] for i in 1:length(p_ini)],t)
    EF2 = (x-x0)*x0_s + (p-p0)*p0_s - Δs

    return [EF1 ; EF2]
end

function Equilibrium_Function(f!::Function,x::Vector{Float64},p::Float64,t::Float64,x0::Vector{Float64},p0::Float64,x0_s::Vector{Float64},p0_s::Float64,Δs::Float64)
    
    dx = zeros(length(x))
    f!(dx,x,p,t)
    EF1 = dx
    EF2 = (x-x0) ⋅ x0_s + (p-p0)*p0_s - Δs

    return [EF1 ; EF2]
end

function Equilibrium_Function(f!::Function,x::Vector{Float64},p::Float64,t::Float64,x0::Vector{Float64},p0::Float64,x0_s::Vector{Float64},p0_s::Float64,Δs::Float64,p_ini::Vector{Float64},indice::Int64)
    
    dx = zeros(length(x))
    f!(dx,x,[i == indice ? p : p_ini[i] for i in 1:length(p_ini)],t)
    EF1 = dx
    EF2 = (x-x0) ⋅ x0_s + (p-p0)*p0_s - Δs

    return [EF1 ; EF2]
end

function Equilibrium_Jacobian(f::Function,x::Float64,p::Float64,t::Float64,x0_s::Float64,p0_s::Float64)

    s = Taylor1(2)
    r = Taylor1([0.0,0.0],2)

    fx = differentiate(f(x+s,p+r,t+r))(0.0)
    fp = differentiate(f(x+r,p+s,t+r))(0.0)
    
    J = [fx fp;
         x0_s p0_s]

    return J

end

function Equilibrium_Jacobian(f::Function,x::Float64,p::Float64,t::Float64,x0_s::Float64,p0_s::Float64,p_ini::Float64,indice::Int64)

    s = Taylor1(2)
    r = Taylor1([0.0,0.0],2)

    fx = differentiate(f(x+s,[i == indice ? p : p_ini[i] for i in 1:length(p_ini)] .+ r,t+r))(0.0)
    fp = differentiate(f(x+r,p+s,t+r))(0.0)
    
    J = [fx fp;
         x0_s p0_s]

    return J

end

function Equilibrium_Jacobian(f!::Function,x::Vector{Float64},p::Float64,t::Float64,x0_s::Vector{Float64},p0_s::Float64)

    n = length(x)

    s = Taylor1(2)
    r = Taylor1([0.0,0.0],2)

    Jx = zeros(n,n)

    for i in 1:n
        for j in 1:n
            dx = [r for i in 1:n]
            f!(dx,x+[j == k ? s : r for k in 1:n],p+r,t+r)
            Jx[i,j] = differentiate(dx[i])(0.0)
        end
    end

    Jp = zeros(n)

    for i in 1:n
        dx = [r for i in 1:n]
        f!(dx,x.+r,p+s,t+r)
        Jp[i] = differentiate(dx[i])(0.0)
    end
    
    J = [Jx Jp;
         transpose(x0_s) p0_s]

    return J

end

function Equilibrium_Jacobian(f!::Function,x::Vector{Float64},p::Float64,t::Float64,x0_s::Vector{Float64},p0_s::Float64,p_ini::Vector{Float64},indice::Int64)

    n = length(x)

    s = Taylor1(2)
    r = Taylor1([0.0,0.0],2)

    Jx = zeros(n,n)

    for i in 1:n
        for j in 1:n
            dx = [r for i in 1:n]
            f!(dx,x+[j == k ? s : r for k in 1:n],[i == indice ? p : p_ini[i] for i in 1:length(p_ini)] .+ r,t+r)
            Jx[i,j] = differentiate(dx[i])(0.0)
        end
    end

    Jp = zeros(n)

    for i in 1:n
        dx = [r for i in 1:n]
        f!(dx,x.+r,[i == indice ? p + s : p_ini[i] + r for i in 1:length(p_ini)] .+ r,t+r)
        Jp[i] = differentiate(dx[i])(0.0)
    end
    
    J = [Jx Jp;
         transpose(x0_s) p0_s]

    return J

end
