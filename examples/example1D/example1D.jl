include("../../src/NumDiffEq.jl");

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
ramas_de_equilibrio = Equilibrium.(g, x_ini, p_ini, t, Δs, p_fin);
LBP = Limit_Points.(g,ramas_de_equilibrio,t,p_fin)
estabilidad = Stability_intervals.(g,t,ramas_de_equilibrio)
end;

#-

LBP = [[(LBP[i][j][1],LBP[i][j][2]) for j in 1:length(LBP[i])] for i in 1:length(LBP)];

#-

begin
    plot(title = "Puntos de equilibrio \n Tiempo de ejecución = $(tiempo) s", ylabel = "u(λ)", xlabel = "λ",legendfont=font(5))
    plot!(estabilidad[1][1], label = "Puntos estables", color = "blue", linestyle = :solid)
    plot!(estabilidad[1][2], label = "Puntos inestables", color = "blue", linestyle = :dash)
    scatter!(LBP[1], label = "Bifurcaciones LP", color = "white")
    for i in 2:length(estabilidad)
        plot!(estabilidad[i][1], label = "", color = "blue", linestyle = :solid)
        plot!(estabilidad[i][2], label = "", color = "blue", linestyle = :dash)
        scatter!(LBP[i], label = "", color = "white")
    end
    savefig("Ejemplo 1D.png")
end;

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