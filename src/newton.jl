# # Método de Newton
# Este archivo contiene las funciones correspondientes la método de Newton
# calculado de forma vectorial en todos los casos.

#-

"""

    Newton(f::Function,x0::Float64,p0::Float64)

Devuelve la raíz `x::Float64` y `p::Float64` más cercana  a los puntos `x0` y `p0` 
de la función `f` relacionada a la ecuación diferencial.
"""
function Newton(f::Function,x0::Float64,p0::Float64)
    x = x0
    p = p0
    i = 1
    while i <= 30 && abs(f(x,p)) > 1.e-16
        df = Gradient(f,x,p)
        x, p = [x,p] - (f(x,p)/norm(df)^2) * df
        i+=1
    end
    return x, p
end

#-

"""

    Newton(f::Function,x0::Float64,p0::Vector{Float64},indice::Int64)

Devuelve la raíz `x::Float64` y `p::Vector{Float64}` más cercana  a los puntos `x0` y `p0`
de la función `f` relacionada a la ecuación diferencial y variando solamente `p0[inidice]`.
"""
function Newton(f::Function,x0::Float64,p0::Vector{Float64},indice::Int64)
    x = x0
    p = p0
    i = 1
    while i <= 30 && abs(f(x,p)) > 1.e-16
        df = Gradient(f,x,p,indice)
        dif_newton = [x; p] - (f(x,p)/norm(df)^2) * df
        x = dif_newton[1]
        p = dif_newton[2:end]
        i+=1
    end
    return x, p
end

#-

"""

    Newton(f!::Function,x0::Vector{Float64},p0::Float64)

Devuelve la raíz `x::Vector{Float64}` y `p::Float64` más cercana  a los puntos `x0` y `p0` 
del sistema de ecuaciones diferenciales asociado a la función `f!`.
"""
function Newton(f!::Function,x0::Vector{Float64},p0::Float64)
    x = x0
    p = p0
    dx = zeros(length(x))
    f!(dx,x,p)
    i = 1
    while i <= 30 && norm(dx) > 1.e-16
        df = Gradient(f!,x,p)
        dif_newton = [x;p] - ((dx ⋅ dx)/norm(df)^2) * df
        x = dif_newton[1:end-1]
        p = dif_newton[end]
        i+=1
        f!(dx,x,p)
    end
    return x, p
end

#-

"""

    Newton(f!::Function,x0::Vector{Float64},p0::Vector{Float64},indice::Int64)

Devuelve la raíz `x::Vector{Float64}` y `p::Vector{Float64}` más cercana  a los puntos `x0` y `p0`
del sistema de ecuaciones diferenciales asociado a la función `f!` y variando solamente `p0[inidice]`.
"""
function Newton(f!::Function,x0::Vector{Float64},p0::Vector{Float64},indice::Int64)
    x = x0
    p = p0
    dx = zeros(length(x))
    f!(dx,x,p)
    i = 1
    while i <= 30 && norm(dx) > 1.e-16
        df = Gradient(f!,x,p,indice)
        dif_newton = [x;p] - ((dx ⋅ dx)/norm(df)^2) * df
        x = dif_newton[1:length(x)]
        p = dif_newton[length(x)+1:end]
        i+=1
        f!(dx,x,p)
    end
    return x, p
end

#-



