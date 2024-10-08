include("../../src/NumDiffEq.jl");

#-

using Plots;

#-

g(u,λ,t) = (u^2-1)*(u^2-4) + λ*u^2*exp(u/10);

#-

begin
    x_ini = [2.0,-1.0]
    p_ini = 0.0
    Δs = 0.001
    t = 0.0
    N = 5000
    indice = 1
end;

#-


tiempo = @elapsed begin
    EB = Equilibrium.(g, x_ini, p_ini, t, Δs, N; bifurcations = true);
end;

#-

begin
    plot(title = "Puntos de equilibrio \n Tiempo de ejecución = $(tiempo) s", ylabel = "u(λ)", xlabel = "λ")

    p, x, estable, inestable, pf, xf = EB[1]

    plot!(p[estable], x[estable], linestyle = :solid, label = "Estable", color = "blue")
    plot!(p[inestable], x[inestable], linestyle = :dash, label = "Inestable", color = "blue")
    scatter!(pf, xf, label = "Bifurcacion LP", color = "white")

    p, x, estable, inestable, pf, xf = EB[2]

    plot!(p[estable], x[estable], linestyle = :solid, label = "", color = "blue")
    plot!(p[inestable], x[inestable], linestyle = :dash, label = "", color = "blue")
    scatter!(pf, xf, label = "", color = "white")

    xlims!(0.0,1.5)

    savefig("Ejemplo 1D.png")
end;

p

#-

begin
    precisiones = []

    #=for i in 1:length(x_ini)
        precisioni = plot(title = "λ_ini = $(p_ini), x_ini = $(x_ini[i])", ylabel = "|g(u(λ),λ)| / 1×10⁻¹⁶", titlefontsize = 7, xlabel = "λ", guidefontsize = 8)
        plot!(ramas_de_equilibrio[i][1],abs.(g.(ramas_de_equilibrio[i][2],ramas_de_equilibrio[i][1],t)) ,label = "")
        push!(precisiones,precisioni)
    end=#

    #plot(precisiones..., layout = (1,2),plot_title="Precisión de las ramas de soluciones",plot_titlefontsize = 10)

    #savefig("Precision.png")
end;

#-