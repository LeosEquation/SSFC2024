

function coefs_taylor(f,t::Taylor1,u::Taylor1,p)

    u_old = copy(u(t))
    u_new = u(t) + integrate(f(u_old,p,t))   
    i = 0
    n = length(u)

    while u_new != u_old

        i += 1
        u_old = u_new
        u_new = u(t) + integrate(f(u_old,p,t))

        if i == (n + 5)
            
            throw(ArgumentError("La iteración de Picard excedió los límites esperados. Verifica tu ecuación diferencial."))
        
        end

    end

    return u_new

end

function coefs_taylor!(f!,t::Taylor1,du::Vector{Taylor1{Float64}},u::Vector{Taylor1{Float64}},p)
    
    u_old = copy(u(t)) 
    f!(du,u_old,p,t)
    u_new = u(t) .+ integrate.(du)
    i = 0
    n = maximum(length.(u))
    m = length(u)

    while u_new != u_old
        
        f!(du,u_new,p,t)
        u_old .= u_new
        u_new .= u(t) .+ integrate.(du)
        i+=1

        if i == (n*m + 5)
            
            throw(ArgumentError("La iteración de Picard excedió los límites esperados. Verifica tu ecuación diferencial."))
        
        end

    end

    return u_new

end



function paso_integracion(u::Taylor1,ϵ::Float64)

    N = (length(u)-1):-1:0
    h = [0.0,0.0]
    j = 0

    for i in N

        if u[i] != 0 && j < 2

            j += 1
            h[j] = (ϵ /abs(u[i]))^(1/i)

        end

        if j == 2

            break

        end

    end

    return 0.5*minimum(h)
end

function paso_integracion(u::Vector{Taylor1{Float64}},ϵ::Float64)

    H = [[0.0,0.0] for i in 1:length(u)]

    for j in 1:length(u)

        k = 0
        N = (length(u[j])-1):-1:0
        
        for i in N

            if u[j][i] != 0 && k < 2

                k += 1
                H[j][k] = (ϵ/abs(u[j][i]))^(1/i)
                
            end

            if k == 2

                break

            end

        end

    end

    return 0.5*minimum(minimum.(H))
    
end

function paso_taylor(f,t::Taylor1,u::Taylor1,p,ϵ::Float64)

    u_sol = coefs_taylor(f,t,u,p)

    return u_sol, paso_integracion(u_sol,ϵ)

end

function paso_taylor(f!,t::Taylor1,du::Vector{Taylor1{Float64}},u::Vector{Taylor1{Float64}},p,ϵ::Float64)

    u_sol = coefs_taylor!(f!,t,du,copy(u),p)

    return u_sol, paso_integracion(u_sol,ϵ)
end

"

# Caso escalar

    t, x, error = integracion_taylor(f, x0::Float64, t_ini::Float64, t_fin::Float64, orden::Int64, ϵ::Float64, p; Nt = 1000, longitud = false)

Este método integra en el caso escalar ecuaciones diferenciales ordinarias de primer orden (ODE's).

`f(x,p,t)` es la función que define la ecuación de la forma `dx/dt = f(x,p,t)`.

`x0::Float64` es la condición inicial de `x` al tiempo `t_ini`, es decir, `x(t = t_ini) = x0`.

`t_ini::Float64` es el tiempo inicial de integración, no es necesariamente menor a `t_fin`.

`t_fin::Float64` es el tiempo final de integración, no es necesariamente mayor a `t_ini`.

`orden::Int64` es el orden máximo del polinomio de Taylor de `x(t)`. Este es siempre positivo. 

`ϵ::Float64` es la tolerancia de error en `x(t)`.

`p` son los parámetros necesarios en la ecuación diferencial `f`. Puede ser un número, un vector si son varios parámetros o puede ser `nothing` si `f` no lo necesita, esto será a elección del usuario.

`t::Vector{Float64}` es un vector de tiempo con valores en el intervalo `[t_ini,t_fin]` incluyendo siempre las fronteras donde la longitud del vector dependerá del `orden` y de la tolerancia `ϵ`. 

`x::Vector{Float64}` es un vector con la solución `x(t)` para todos los valores del vector de tiempo `t` dentro de la tolerancia `ϵ`.

`error::Bool` es una variable Booleana. Si `error = true` indica que el error de la solución entra dentro de la tolerancia `ϵ`, si `error = false` indica el caso contrario.

### Precaución con los parámetros opcionales `Nt::Int64` y `longitud::Bool`


Si `longitud = true` el intervalo `[t_ini,t_fin]` se partirá forsoza y homogeneamente en `Nt` partes sin importar los valores de `ϵ` y `orden`.

Si la tolerancia `ϵ` es muy pequeña y el `orden` es pequeño, el intervalo `[t_ini,t_fin]` se partirá homogeneamente en `Nt` partes sin importar el valor de `longitud` pero la tolerancia ya no será respetada, es decir, el error en `x` será mayor al de la tolerancia `ϵ`.

### Ejemplo de una ecuación diferencial

Para definir la ecuación diferencial `dx/dt = 2.0*x`


    julia> f(x,p,t) = p*x
    f (generic function with 1 method)

    julia> p = 2.0
    2.0


O también puede ser


    julia> f(x,p,t) = 2.0*x
    f (generic function with 1 method)

    julia> p = nothing


Si son varios parámetros, por ejemplo `dx/dt = -3.0 + 2.0*x - 5.0*x^2`, p se puede implementar así


    julia> f(x,p,t) = p[1] + p[2]*x + p[3]*x^2
    f (generic function with 1 method)

    julia> p = [-3.0,2.0,-5.0]
    3-element Vector{Float64}:
     -3.0
      2.0
     -5.0
"
function integracion_taylor(f, x0::Float64, t_ini::Float64, t_fin::Float64, orden::Int64, ϵ::Float64, p; Nt = 1000, longitud = false)

    u = Taylor1(x0,orden)
    t = Taylor1([t_ini,1],orden)
    u, h = paso_taylor(f,t,copy(u),p,ϵ)
    t_min = minimum([t_ini,t_fin])
    t_max = maximum([t_ini,t_fin])
    Δt = t_max - t_min
    tt = t_ini
    T = [t_ini]
    X = [x0]

    Δh = Δt/Nt

    if h < Δh || longitud
        h = Δh
    end

    hh = (-1)^(t_ini > t_fin) * h

    while t_min <= tt + hh <= t_max

        tt += hh
        u = Taylor1(copy(u(hh)),orden)
        t = Taylor1([tt,1],orden)
        u, h = paso_taylor(f,t,copy(u),p,ϵ)
        push!(T,tt)
        push!(X,u(0))

        if h < Δh || longitud
             h = Δh
        end

        hh = (-1)^(t_ini > t_fin) * h

    end

    u = Taylor1(copy(u(t_fin - tt)),orden)
    t = Taylor1([t_fin,1],orden)
    u = coefs_taylor(f,t,copy(u),p)
    push!(T,t_fin)
    push!(X,u(0))

    return T, X

end

"

# Caso vectorial


    t, x, error = integracion_taylor(f!, x0::Vector{Float64}, t_ini::Float64, t_fin::Float64, orden::Int64, ϵ::Float64, p; N = 1000, longitud = false)

Este método integra en el caso vectorial ecuaciones diferenciales ordinarias de primer orden (ODE's).

`f!(dx,x,p,t)` es la función que define la ecuación vectorial. Esta función cambiará la forma de `dx = [dx_1/t,dx_2/t,dx_3/t,...]` que se define automáticamente dentro de `integracion_taylor`.

`x0::Vector{Float64}` es la condición inicial del vector `x` al tiempo `t_ini`, es decir, `x(t = t_ini) = x0 = (x10,x20,x30,...)`.

`x::Vector{Vector{Float64}}` es un vector con los vectores solución de `x(t)` para todos los valores del vector de tiempo `t` dentro de la tolerancia `ϵ`.

### Ejemplo de una ecuación diferencial vectorial

Para definir la ecuación diferencial vectorial `dx/dt = (dx1/dt,dx2/dt) = (x2,-x1)`


    julia> function f!(dx,x,p,t)
                dx[1] = x[2]
                dx[2] = -x[1]
           end
    f (generic function with 1 method)

"
function integracion_taylor(f!, x0::Vector{Float64}, t_ini::Float64, t_fin::Float64, orden::Int64, ϵ::Float64, p; Nt = 1000)

    u = Taylor1.(x0,orden)  
    t = Taylor1([t_ini,1],orden) 
    du = Taylor1.(zero(x0),orden)
    uu, h = paso_taylor(f!,t,du,copy(u),p,ϵ)
    u .= uu
    t_min = minimum([t_ini,t_fin])
    t_max = maximum([t_ini,t_fin])
    Δt = t_max - t_min
    signo = sign(t_fin-t_ini)
    tt = t_ini
    Δh = Δt/Nt

    if t_ini < t_fin
        T = t_ini:Δh:t_fin
    else
        T = t_ini:-Δh:t_fin
    end
    #T = [tt]
    X = [x0]

    #Δh = Δt/Nt

    #if h < Δh || longitud
        #h = Δh
    #end

    #hh = (-1)^(t_ini > t_fin) * h
    hh = signo * h
    i = 1
    j = 2

    while t_min <= tt + hh <= t_max

        tt += hh
        u .= Taylor1.(copy(u(hh)),orden)
        t = Taylor1([tt,1],orden)
        uu, h = paso_taylor(f!,t,du,copy(u),p,ϵ)
        u .= uu
        #push!(T,tt)
        while (tt - hh) <= T[j] < tt && j < length(T)
            push!(X,u(T[j]-tt))
            j += 1
        end

        #for i in 1:(length(T)-1)
            #if tt <= T[i] < (tt + hh)
                #push!(X,u(T[j]-tt))
            #end
        #end

        #if h < Δh || longitud
            #h = Δh
        #end

        #hh = (-1)^(t_ini > t_fin) * h

        if h < Δh
            h = Δh
        end
        
        hh = signo * h

    end

    u .= Taylor1.(copy(u(t_fin - tt)),orden)
    t = Taylor1([t_fin,1],orden)
    u .= coefs_taylor!(f!,t,du,copy(u),p)
    #push!(T,t_fin)
    for k in (j):length(T)
        push!(X,u(T[k]-t_fin))
    end

    return T, X

end
