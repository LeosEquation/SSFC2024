# # Bifurcaciones
# Este archivo contiene las funciones correspondientes para encontrar
# puntos de bifurcación en un sistema de ecuaciones diferenciales.

#-

function Limit_Points(f::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Float64}},t::Float64,p_fin::Float64; ω = 1.0, test = false)
    
    initial_values = LP_initial_values(f,equilibrium_branch,t)

    p_min = minimum([equilibrium_branch[1][1],p_fin])
    p_max = maximum([equilibrium_branch[1][1],p_fin])

    LPs = Vector{Float64}[]

    for initial_value in initial_values
        p0, x0, v0, λ0 = initial_value

        x = x0
        p = p0
        v = v0

        i = 1

        while i <= 30 && norm(LP_Function(f,x,p,v,t; ω = ω)) > 1.e-16 && norm(λ0) > 1.e-16

            if test
                println("Paso $(i): λ = $(p)")
                println("           x = $(x)")
                println("           norm(LP_Function) = $(norm(LP_Function(f,x,p,v,t)))")
                println("\n")
                show(stdout, "text/plain",  LP_Function(f,x,p,v,t))
                println("\n")
            end

            ω = v

            X = [x;p;v] - inv(LP_Jacobian(f,x,p,v,t; ω = ω))*LP_Function(f,x,p,v,t; ω = ω)
            x = X[1]
            p = X[2]
            v = X[3]
            i += 1
        end
        
        if p_min <= p <= p_max
            push!(LPs,[p,x])
        end

    end

        return LPs
end

#-

function Limit_Points(f::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Float64}},t::Float64,p_ini::Vector{Float64},indice::Int64,p_fin::Float64; ω = 1.0, test = false)
    
    initial_values = LP_initial_values(f,equilibrium_branch,t,p_ini,indice)

    p_min = minimum([equilibrium_branch[1][1],p_fin])
    p_max = maximum([equilibrium_branch[1][1],p_fin])

    LPs = Vector{Float64}[]

    for initial_value in initial_values
        p0, x0, v0, λ0 = initial_value

        x = x0
        p = p0
        v = v0

        i = 1

        while i <= 30 && norm(LP_Function(f!,x,p,v,t,p_ini,indice; ω = ω)) > 1.e-16 && norm(λ0) > 1.e-16

            if test
                println("Paso $(i): λ = $(p)")
                println("           x = $(x)")
                println("           v = $(v)")
                println("           norm(LP_Function) = $(norm(LP_Function(f!,x,p,v,t,p_ini,indice; ω = ω)))")
                println("\n")
                show(stdout, "text/plain",  LP_Function(f!,x,p,v,t,p_ini,indice; ω = ω))
                println("\n")
            end

            ω = v

            X = [x;p;v] - inv(LP_Jacobian(f!,x,p,v,t,p_ini,indice; ω = ω))*LP_Function(f!,x,p,v,t,p_ini,indice; ω = ω)
            x = X[1]
            p = X[2]
            v = X[3]
            i += 1
        end
        
        if p_min <= p <= p_max
            push!(LPs,[p,x])
        end

    end

        return LPs
end
#-

function Limit_Points(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},t::Float64,p_fin::Float64; ω = ones(length(equilibrium_branch[2][1])), test = false)
    
    initial_values = LP_initial_values(f!,equilibrium_branch,t)

    p_min = minimum([equilibrium_branch[1][1],p_fin])
    p_max = maximum([equilibrium_branch[1][1],p_fin])

    LPs = Vector{Union{Float64,Vector{Float64}}}[]

    for initial_value in initial_values
        p0, x0, v0, λ0 = initial_value

        x = x0
        p = p0
        v = v0

        i = 1

        while i <= 30 && norm(LP_Function(f!,x,p,v,t; ω = ω)) > 1.e-16 && norm(λ0) > 1.e-16

            if test
                println("Paso $(i): λ = $(p)")
                println("           x = $(x)")
                println("           norm(LP_Function) = $(norm(LP_Function(f!,x,p,v,t)))")
                println("\n")
                show(stdout, "text/plain",  LP_Function(f!,x,p,v,t))
                println("\n")
            end

            ω = v

            X = [x;p;v] - inv(LP_Jacobian(f!,x,p,v,t; ω = ω))*LP_Function(f!,x,p,v,t; ω = ω)
            x = X[1:length(x)]
            p = X[length(x)+1]
            v = X[length(x)+2:end]
            i += 1
        end
        
        if p_min <= p <= p_max
            push!(LPs,[p,x])
        end

    end

        return LPs
end

#-

function Limit_Points(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},t::Float64,p_ini::Vector{Float64},indice::Int64,p_fin::Float64; ω = ones(length(equilibrium_branch[2][1])), test = false)
    
    initial_values = LP_initial_values(f!,equilibrium_branch,t,p_ini,indice)

    p_min = minimum([equilibrium_branch[1][1],p_fin])
    p_max = maximum([equilibrium_branch[1][1],p_fin])

    LPs = Vector{Union{Float64,Vector{Float64}}}[]

    for initial_value in initial_values
        p0, x0, v0, λ0 = initial_value

        x = x0
        p = p0
        v = v0

        i = 1

        while i <= 30 && norm(LP_Function(f!,x,p,v,t,p_ini,indice; ω = ω)) > 1.e-16 && norm(λ0) > 1.e-16

            if test
                println("Paso $(i): λ = $(p)")
                println("           x = $(x)")
                println("           v = $(v)")
                println("           norm(LP_Function) = $(norm(LP_Function(f!,x,p,v,t,p_ini,indice; ω = ω)))")
                println("\n")
                show(stdout, "text/plain",  LP_Function(f!,x,p,v,t,p_ini,indice; ω = ω))
                println("\n")
            end

            ω = v

            X = [x;p;v] - inv(LP_Jacobian(f!,x,p,v,t,p_ini,indice; ω = ω))*LP_Function(f!,x,p,v,t,p_ini,indice; ω = ω)
            x = X[1:length(x)]
            p = X[length(x)+1]
            v = X[length(x)+2:end]
            i += 1
        end
        
        if p_min <= p <= p_max
            push!(LPs,[p,x])
        end

    end

        return LPs
end

#-

function Hopf_Points(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},t::Float64,p_fin::Float64; test = false)
    
    initial_values = Hopf_initial_values(f!,equilibrium_branch,t)

    p_min = minimum([equilibrium_branch[1][1],p_fin])
    p_max = maximum([equilibrium_branch[1][1],p_fin])

    HPs = Vector{Union{Float64,Vector{Float64}}}[]

    for initial_value in initial_values
        p0, x0, κ0, v0, w0 = initial_value

        x = x0
        p = p0
        v = v0
        κ = κ0
        w = w0
        
        i = 1

        while i <= 30 && norm(Hopf_Function(f!,x,p,v,κ,t,w)) > 1.e-16

            if test
                println("Paso $(i): λ = $(p)")
                println("           x = $(x)")
                println("           norm(Hopf_Function) = $(norm(Hopf_Function(f!,x,p,v,κ,t,w)))")
                println("\n")
                show(stdout, "text/plain",  Hopf_Function(f!,x,p,v,κ,t,w))
                println("\n")
            end

            X = [x;p;v;κ] - inv(Hopf_Jacobian(f!,x,p,v,κ,t,w))*Hopf_Function(f!,x,p,v,κ,t,w)
            x = X[1:length(x)]
            p = X[length(x)+1]
            v = X[length(x)+2:2*length(x)+1]
            κ = X[end]

            i += 1
        end

        if p_min <= p <= p_max
            push!(HPs,[p,x])
        end

    end

        return HPs
end

#-

function Hopf_Points(f!::Function,equilibrium_branch::Tuple{Vector{Float64}, Vector{Vector{Float64}}},t::Float64,p_ini::Vector{Float64},indice::Int64,p_fin::Float64; test = false)
    
    initial_values = Hopf_initial_values(f!,equilibrium_branch,t,p_ini,indice)

    p_min = minimum([equilibrium_branch[1][1],p_fin])
    p_max = maximum([equilibrium_branch[1][1],p_fin])

    HPs = Vector{Union{Float64,Vector{Float64}}}[]

    for initial_value in initial_values
        p0, x0, κ0, v0, w0 = initial_value

        x = x0
        p = p0
        v = v0
        κ = κ0
        w = w0
        
        i = 1

        while i <= 30 && norm(Hopf_Function(f!,x,p,v,κ,t,w,p_ini,indice)) > 1.e-16

            if test
                println("Paso $(i): λ = $(p)")
                println("           x = $(x)")
                println("           norm(Hopf_Function) = $(norm(Hopf_Function(f!,x,p,v,κ,t,w,p_ini,indice)))")
                println("\n")
                show(stdout, "text/plain",  Hopf_Function(f!,x,p,v,κ,t,w,p_ini,indice))
                println("\n")
            end

            X = [x;p;v;κ] - inv(Hopf_Jacobian(f!,x,p,v,κ,t,w,p_ini,indice))*Hopf_Function(f!,x,p,v,κ,t,w,p_ini,indice)
            x = X[1:length(x)]
            p = X[length(x)+1]
            v = X[length(x)+2:2*length(x)+1]
            κ = X[end]

            i += 1
        end
        
        if p_min <= p <= p_max
            push!(HPs,[p,x])
        end

    end

        return HPs
end

#-


# ### Referencias

#- 

# - http://www.maths.liv.ac.uk/~bnvasiev/Past%20students/Caitlin_399.pdf
# 
# - Doedel, E.J. (2007). Lecture Notes on Numerical Analysis of Nonlinear Equations.
#   In: Krauskopf, B., Osinga, H.M., Galán-Vioque, J. (eds) Numerical Continuation 
#   Methods for Dynamical Systems. Understanding Complex Systems. Springer, Dordrecht. 
#   https://doi.org/10.1007/978-1-4020-6356-5_1 