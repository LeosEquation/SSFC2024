include("../../src/NumDiffEq.jl");

using .NumDiffEq
using Plots, TaylorSeries, LinearAlgebra

#-

function f!(du,u,λ,t)
    du[1] = 3*u[1]*(1 - u[1]) - u[1]*u[2] - λ*(1-exp(-5*u[1]))
    du[2] = -u[2] + 3*u[1]*u[2]
end;

#-

begin
    x_ini = [[0.0,0.0],[1.0,0.0],[1/3,2.0]]
    p_ini = 0.0
    Δs= 0.001
    p_fin = 0.9
    t = 0.0
end;

#-

tiempo = @elapsed begin
    ramas_de_equilibrio = Equilibrium.(f!, x_ini, p_ini, t, Δs, p_fin)
    LBP = Limit_Points.(f!,ramas_de_equilibrio,Δs,t)
    estabilidad = Stability_intervals.(f!,t,ramas_de_equilibrio)
end;

#-

begin
    LBP_2D = [[(LBP[i][j][1],LBP[i][j][2][1]) for j in 1:length(LBP[i])] for i in 1:length(LBP)]
    estable = [[(estabilidad[i][1][j][1],estabilidad[i][1][j][2][1]) for j in 1:length(estabilidad[i][1])] for i in 1:length(estabilidad)]
    inestable = [[(estabilidad[i][2][j][1],estabilidad[i][2][j][2][1]) for j in 1:length(estabilidad[i][2])] for i in 1:length(estabilidad)]
end;

#-

begin
    plot(title = "Predator-Prey Model \n Tiempo de ejecución = $(tiempo) s", ylabel = "Fish", xlabel = "Quota")
    plot!(estable[1], label = "Estable", linestyle = :solid, color = "blue")
    plot!(inestable[1], label = "Inestable", linestyle = :dash, color = "blue")
    scatter!(LBP_2D[1], label = "Bifurcación LP", color = "white")
    ylims!(-0.1,1.0)
    xlims!(0.0,0.9)

    for i in 1:length(ramas_de_equilibrio)
    plot!(estable[i], label = "", linestyle = :solid, color = "blue")
    plot!(inestable[i], label = "", linestyle = :dash, color = "blue")
    scatter!(LBP_2D[i], label = "", color = "white")
    end

    savefig("Predator-Prey fish.png")
end;

#-

begin
    LBP_3D = [[(LBP[i][j][2][2],LBP[i][j][1],LBP[i][j][2][1]) for j in 1:length(LBP[i])] for i in 1:length(LBP)]

    plot(title = "Predator-Prey Model \n Tiempo de ejecución = $(tiempo) s", ylabel = "Quota", xlabel = "Sharks", zlabel = "Fish", camera = (60, 25))

    for i in 1:3
        x = zeros(length(ramas_de_equilibrio[i][1]))

        y_estable = zeros(length(ramas_de_equilibrio[i][1]))
        z_estable = zeros(length(ramas_de_equilibrio[i][1]))

        y_inestable = zeros(length(ramas_de_equilibrio[i][1]))
        z_inestable = zeros(length(ramas_de_equilibrio[i][1]))

        x = ramas_de_equilibrio[i][1]

        for j in 1:length(estabilidad[i][1])
            y_estable[j] = estabilidad[i][1][j][2][2]
            z_estable[j] = estabilidad[i][1][j][2][1]
        end

        for j in 1:length(estabilidad[i][2])
            y_inestable[j] = estabilidad[i][2][j][2][2]
            z_inestable[j] = estabilidad[i][2][j][2][1]
        end

        if i == 1
            plot!(y_estable,x,z_estable,xflip = true, color = "blue",linestyle = :solid, label = "Estable")
            plot!(y_inestable,x,z_inestable,xflip = true, color = "blue",linestyle = :dash, label = "Inestable")
            scatter!(LBP_3D[i], label = "Bifurcación LP", color = "white")
        else
            plot!(y_estable,x,z_estable,xflip = true, color = "blue",linestyle = :solid, label = "")
            plot!(y_inestable,x,z_inestable,xflip = true, color = "blue",linestyle = :dash, label = "")
            scatter!(LBP_3D[i], label = "", color = "white")
        end
    end


    ylims!(-0.1,1.0)
    xlims!(-0.5,2.0)
    zlims!(-0.25,1.0)
    savefig("Predator-Prey Model.png")
end

#-

begin
    for i in 1:length(x_ini)
        normas = []
        for j in 1:length(solfams[i][1])
            dx = zeros(length(x_ini))
            f!(dx,solfams[i][2][j],solfams[i][1][j],t)
            push!(normas,norm(dx))
        end
        if i == 2
            precision = plot(title = "Precisión de la rama", ylabel = "||F(u(λ),λ)||", xlabel = "λ")
            plot!(solfams[i][1],normas ,label = "λ_ini = $(p_ini), x_ini = $(x_ini[i])")
        else
            precision = plot(title = "Precisión de la rama", ylabel = "||F(u(λ),λ)|| / 1×10⁻¹⁶", xlabel = "λ")
            plot!(solfams[i][1],normas .* 1.e16 ,label = "λ_ini = $(p_ini), x_ini = $(x_ini[i])")
        end
        savefig("PrecisiónRama$(i).png")
    end
end

#-

function ν_initial(f!,t,x0,p0)

    n = length(x0)
    ceros = zeros(n+1)

    S = set_variables("s", numvars = n+1, order = 3)

    x = x0 + S[1:end-1]
    p = p0 + S[end]

    dx = x

    f!(dx,x,p,t)

    Fx = [derivative(dx[i],j)(ceros) for i in 1:n, j in 1:n]

    ν = [1.0,1.0]
    ω = imag(eigvals(Fx)[1])

    for i in 1:100

        μ = (Fx^2 + ω^2*I(n))\ν

        ν = μ/norm(μ)

    end

    return ν

end

prueba = Equilibrium(f!, x_ini[3], p_ini, t, Δs, p_fin; N = 10000)

p0, x0, ν0, λ0, w0 = Hopf_initial_values(f!,prueba,t)[1]

x = x0
p = p0
κ = imag(λ0[1])^2
ν = real(ν0[:,1])
w = w0[:,1]

ν = ν_initial(f!,t,x,p)
w = [ν[2],-ν[1]]

ν ⋅ w

Hopf_function(f!,x,p,ν,ω^2,t,w)

i = 1

while i <= 30 && norm(Hopf_function(f!,x,p,ν,κ,t,w)) > 1.e-16
    println("Paso $(i): λ = $(p)")
    println("         x = $(x)")
    println("         ν = $(ν)")
    println("         κ = $(κ)")
    println("         w = $(w)")
    println("         norm(Hopf_function) = $(norm(Hopf_function(f!,x,p,ν,κ,t,w)))")
    println("\n")
    show(stdout, "text/plain",  Hopf_function(f!,x,p,ν,κ,t,w))
    println("\n")

    X = [x;p;ν;κ] - inv(Hopf_Jacobian(f!,x,p,ν,κ,t,w))*Hopf_function(f!,x,p,ν,κ,t,w)
    x = X[1:length(x)]
    p = X[length(x)+1]
    ν = X[length(x)+2:2*length(x)+1]
    κ = X[end]

    i += 1
end

norm(Hopf_function(f!,x,p,ν,κ,t,w))

x
p
ν
κ

