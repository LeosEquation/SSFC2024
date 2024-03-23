# # Pseudo Arc Length Continuation
# Este archivo contiene las funciones que calculan los pasos
# en el método Pseudo Arc Length Continuation

#-

"""

    Derivative_arclength(f::Function,x_ini::Float64,p_ini::Float64,p_fin::Float64)

Devuelve las derivadas `x_s::Float64, p_s::Float64` de las variables `x` y `p` respecto 
a la longitud de arco `s` en `x_ini` y `p_ini`.
"""
function Derivative_arclength(f::Function,x_ini::Float64,p_ini::Float64,p_fin::Float64)

    x_s = 1/sqrt(1/(Derivative_IFT(f,x_ini,p_ini)^2) + 1)
    p_s = sign(p_fin - p_ini)/sqrt(Derivative_IFT(f,x_ini,p_ini)^2 + 1)

    return x_s, p_s

end

#-

"""

    step(f::Function,x::Float64, p::Float64, x_s::Float64,p_s::Float64)

Devuelve una aproximación numérica de las siguientes derivadas `x_s_new::Float64, p_s_new::Float64`
respecto a la longitud de arco `s` evaluadas en `x` y `p`.
"""
function step(f::Function,x::Float64, p::Float64, x_s::Float64,p_s::Float64)
    t = Taylor1(1)
    
    f_x = derivative(f(x+t,p))(0.0)
    f_p = derivative(f(x,p+t))(0.0)

    A = [f_x f_p;
        x_s p_s]

    b = transpose([0 1])

    x_s_new, p_s_new = A\b

    return x_s_new, p_s_new
end

#-

"""

    Derivative_arclength(f::Function,x_ini::Float64,p_ini::Vector{Float64},p_fin::Float64,indice::Int64)

Devuelve las derivadas `x_s::Float64, p_s::Float64` de las variables `x` y `p[indice]` respecto 
a la longitud de arco `s` en `x_ini` y `p_ini[indice]`.
"""
function Derivative_arclength(f::Function,x_ini::Float64,p_ini::Vector{Float64},p_fin::Float64, indice::Int64)
    x_s = 1/sqrt(1/(Derivative_IFT(f,x_ini,p_ini,indice)^2) + 1)
    p_s = sign(p_fin - p_ini[indice])/sqrt(Derivative_IFT(f,x_ini,p_ini,indice)^2 + 1)

    return x_s, p_s

end

#-

"""

    step(f::Function,x::Float64, p::Float64, x_s::Float64,p_s::Vector{Float64},indice::Int64)

Devuelve una aproximación numérica de las siguientes derivadas `x_s_new::Float64, p_s_new::Float64`
respecto a la longitud de arco `s` evaluadas en `x` y `p[ini]`.
"""
function step(f::Function,x_ini::Float64, p_ini::Vector{Float64}, x_s::Float64,p_s::Float64,indice::Int64)
    t = Taylor1(1)
    T = [i == indice ? t : Taylor1(0) for i in 1:length(p_ini)]

    f_x = derivative(f(x_ini+t,p_ini))(0.0)
    f_p = derivative(f(x_ini,p_ini+T))(0.0)

    A = [f_x f_p;
        x_s p_s]

    b = transpose([0 1])

    x_s_new, p_s_new = A\b

    return x_s_new, p_s_new
end

#-

"""

    Derivative_arclength(f::Function,x_ini::Vector{Float64},p_ini::Float64,p_fin::Float64)

Devuelve las derivadas `x_s::Vector{Float64}, p_s::Float64` de las variables `x` y `p` respecto 
a la longitud de arco `s` en `x_ini` y `p_ini`.
"""
function Derivative_arclength(f!::Function,x_ini::Vector{Float64},p_ini::Float64,p_fin::Float64)

    p_s = sign(p_fin - p_ini)/sqrt(norm(Derivative_IFT(f!,x_ini,p_ini))^2 + 1)
    x_s = abs(p_s)*Derivative_IFT(f!,x_ini,p_ini)
    return x_s, p_s

end

#-

"""

    step(f::Function,x::Vector{Float64}, p::Float64, x_s::Float64,p_s::Float64)

Devuelve una aproximación numérica de las siguientes derivadas `x_s_new::Vector{Float64}, p_s_new::Float64`
respecto a la longitud de arco `s` evaluadas en `x` y `p`.
"""
function step(f!::Function,x::Vector{Float64}, p::Float64, x_s::Vector{Float64},p_s::Float64)
    
    J = Jacobian(f!,x,p)

    A = [J;
        transpose(x_s) p_s]

    m = length(x_s)+1

    b = [i != m ? 0.0 : 1.0 for i in 1:m]

    Solution = A\b

    #println(Solution)

    return Solution[1:end-1], Solution[end]
end

#-

"""

    Derivative_arclength(f::Function,x_ini::Vector{Float64},p_ini::Vector{Float64},p_fin::Float64,indice::Int64)

Devuelve las derivadas `x_s::Vector{Float64}, p_s::Float64` de las variables `x` y `p[indice]` respecto 
a la longitud de arco `s` en `x_ini` y `p_ini[indice]`.
"""
function Derivative_arclength(f!::Function,x_ini::Vector{Float64},p_ini::Vector{Float64},p_fin::Float64, indice::Int64)
    p_s = sign(p_fin - p_ini[indice])/sqrt(norm(Derivative_IFT(f!,x_ini,p_ini,indice))^2 + 1)
    x_s = abs(p_s)*Derivative_IFT(f!,x_ini,p_ini,indice)
    return x_s, p_s

end

#-

"""

    step(f::Function,x::Vector{Float64}, p::Float64, x_s::Vector{Float64},p_s::Vector{Float64},indice::Int64)

Devuelve una aproximación numérica de las siguientes derivadas `x_s_new::Vector{Float64}, p_s_new::Float64`
respecto a la longitud de arco `s` evaluadas en `x` y `p[ini]`.
"""

function step(f!::Function,x::Vector{Float64}, p::Vector{Float64}, x_s::Vector{Float64},p_s::Float64,indice::Int64)

    J = Jacobian(f!,x,p,indice)

    A = [J;
        transpose(x_s) p_s]

    m = length(x_s)+1

    b = [i != m ? 0.0 : 1.0 for i in 1:m]

    Solution = A\b

    #println(Solution)

    return Solution[1:end-1], Solution[end]
end

#-

