# # Pseudo Arc Length Continuation
# Este archivo contiene las funciones que calculan los pasos
# en el método Pseudo Arc Length Continuation

#-

"""

    Derivative_arclength(f::Function,x_ini::Float64,p_ini::Float64,p_fin::Float64)

Devuelve las derivadas `x_s::Float64, p_s::Float64` de las variables `x` y `p` respecto 
a la longitud de arco `s` en `x_ini` y `p_ini`.
"""
function Derivative_arclength(f::Function,x_ini::Float64,p_ini::Float64,t::Float64,p_fin::Float64)
    
    p_s = 1/sqrt(Derivative_IFT(f,x_ini,p_ini,t)^2 + 1)
    x_s = p_s*Derivative_IFT(f,x_ini,p_ini,t)
    

    return x_s, sign(p_fin - p_ini)*p_s

end

#-

"""

    step(f::Function,x::Float64, p::Float64, x_s::Float64,p_s::Float64)

Devuelve una aproximación numérica de las siguientes derivadas `x_s_new::Float64, p_s_new::Float64`
respecto a la longitud de arco `s` evaluadas en `x` y `p`.
"""
function step(f::Function,x::Float64, p::Float64, t::Float64, x_s::Float64,p_s::Float64)
    
    f_x = derivative(f(x + s,p + r, t + r))(0.0)
    f_p = derivative(f(x + r,p + s, t + r))(0.0)

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
function Derivative_arclength(f::Function,x_ini::Float64,p_ini::Vector{Float64},t::Float64,p_fin::Float64, indice::Int64)
    
    p_s = 1/sqrt(Derivative_IFT(f,x_ini,p_ini,t,indice)^2 + 1)
    x_s = p_s * Derivative_IFT(f,x_ini,p_ini,t,indice)
    #x_s = 1/sqrt(1/(Derivative_IFT(f,x_ini,p_ini,t,indice)^2) + 1)
    
    return x_s, sign(p_fin - p_ini[indice])*p_s

end

#-

"""

    step(f::Function,x::Float64, p::Float64, x_s::Float64,p_s::Vector{Float64},indice::Int64)

Devuelve una aproximación numérica de las siguientes derivadas `x_s_new::Float64, p_s_new::Float64`
respecto a la longitud de arco `s` evaluadas en `x` y `p[ini]`.
"""
function step(f::Function,x_ini::Float64, p_ini::Vector{Float64},t::Float64, x_s::Float64,p_s::Float64,indice::Int64)

    f_x = derivative(f(x_ini + s, p_ini .+ r, t + r))(0.0)
    f_p = derivative(f(x_ini .+ r, p_ini + S, t + r))(0.0)

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
function Derivative_arclength(f!::Function,x_ini::Vector{Float64},p_ini::Float64, t::Float64,p_fin::Float64)

    p_s = 1/sqrt(norm(Derivative_IFT(f!,x_ini,p_ini,t))^2 + 1)
    x_s = p_s*Derivative_IFT(f!,x_ini,p_ini,t)
    return x_s, sign(p_fin - p_ini)*p_s

end

#-

"""

    step(f::Function,x::Vector{Float64}, p::Float64, x_s::Float64,p_s::Float64)

Devuelve una aproximación numérica de las siguientes derivadas `x_s_new::Vector{Float64}, p_s_new::Float64`
respecto a la longitud de arco `s` evaluadas en `x` y `p`.
"""
function step(f!::Function,x::Vector{Float64}, p::Float64, t::Float64, x_s::Vector{Float64},p_s::Float64)
    
    J = Jacobian(f!,x,p,t)

    A = [J;
        transpose(x_s) p_s]

    m = length(x_s)+1

    b = [i != m ? 0.0 : 1.0 for i in 1:m]

    Solution = A\b

    return Solution[1:end-1], Solution[end]
end

#-

"""

    Derivative_arclength(f::Function,x_ini::Vector{Float64},p_ini::Vector{Float64},p_fin::Float64,indice::Int64)

Devuelve las derivadas `x_s::Vector{Float64}, p_s::Float64` de las variables `x` y `p[indice]` respecto 
a la longitud de arco `s` en `x_ini` y `p_ini[indice]`.
"""
function Derivative_arclength(f!::Function,x_ini::Vector{Float64},p_ini::Vector{Float64},t::Float64,p_fin::Float64, indice::Int64)
    p_s = 1/sqrt(norm(Derivative_IFT(f!,x_ini,p_ini,t,indice))^2 + 1)
    x_s = p_s*Derivative_IFT(f!,x_ini,p_ini,t,indice)
    return x_s, sign(p_fin - p_ini[indice])*p_s

end

#-

"""

    step(f::Function,x::Vector{Float64}, p::Float64, x_s::Vector{Float64},p_s::Vector{Float64},indice::Int64)

Devuelve una aproximación numérica de las siguientes derivadas `x_s_new::Vector{Float64}, p_s_new::Float64`
respecto a la longitud de arco `s` evaluadas en `x` y `p[ini]`.
"""

function step(f!::Function,x::Vector{Float64}, p::Vector{Float64}, t::Float64, x_s::Vector{Float64},p_s::Float64,indice::Int64)

    J = Jacobian(f!,x,p,t,indice)

    A = [J;
        transpose(x_s) p_s]

    m = length(x_s)+1

    b = [i != m ? 0.0 : 1.0 for i in 1:m]

    Solution = A\b

    return Solution[1:end-1], Solution[end]
end

#-

