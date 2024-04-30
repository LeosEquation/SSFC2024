include("../../src/equilibrium.jl");
include("../../src/bifurcation.jl");
include("../../src/stability.jl");

#-

using Plots;

#-

g(u,λ,t) = (u^2-1)*(u^2-4) + λ*u^2*exp(u/10);

#-

begin
    x_ini = [2.0,-1.0]
    p_ini = 0.0
    Δs= 0.001
    p_fin = 2.0
    orden = 1
    t = 0.0
end;

#-

tiempo = @elapsed begin
solfams = Equilibrium.(g, x_ini, p_ini, t, Δs, p_fin);
bifurcaciones = Bifurcation_point.(g,x_ini,p_ini,t)
estabilidad = Stability_intervals.(g,t,solfams)
end;

#-

begin
    plot(title = "Puntos de equilibrio \n Tiempo de ejecución = $(tiempo) s", ylabel = "u(λ)", xlabel = "λ",legendfont=font(5))
    plot!(estabilidad[1][1], label = "Puntos estables", color = "blue", linestyle = :solid)
    plot!(estabilidad[1][2], label = "Puntos inestables", color = "blue", linestyle = :dash)
    for i in 2:length(estabilidad)
        plot!(estabilidad[i][1], label = "", color = "blue", linestyle = :solid)
        plot!(estabilidad[i][2], label = "", color = "blue", linestyle = :dash)
    end
    scatter!(bifurcaciones[1],color = "red", label = "Bifurcación", markersize = 3)
    scatter!(bifurcaciones[2:end],color = "red", label = "", markersize = 3)
    savefig("Ejemplo 1D.png")
end;

#-

begin
    precisiones = []

    for i in 1:length(x_ini)
        precisioni = plot(title = "λ_ini = $(p_ini), x_ini = $(x_ini[i])", ylabel = "|g(u(λ),λ)| / 1×10⁻¹⁶", titlefontsize = 7, xlabel = "λ", guidefontsize = 8)
        plot!(solfams[i][1],abs.(g.(solfams[i][2],solfams[i][1],t)) .* 1.e16 ,label = "")
        push!(precisiones,precisioni)
    end

    plot(precisiones..., layout = (1,2),plot_title="Precisión de las ramas de soluciones",plot_titlefontsize = 10)

    savefig("Precision.png")
end;

#-

A = [rand() for i in 1:3,j in 1:3]

det(A)

L, U, P = factorize(A)


L
U
P

A

P

[i == j ? 1 : 0 for i in 1:length(P), j in P]*L*U

-det(U)

det([i == j ? 1 : 0 for i in 1:length(P), j in P])

det(L)

det(A)