include("../../src/solution_family.jl");
using Plots, Statistics;

#-

g(u,λ,t) = (u^2-1)*(u^2-4) + λ*u^2*exp(u/10);

#-

x_ini = [2.0,1.0,-1.0,-2.0];
p_ini = 0.0;
Δs= 0.001;
p_fin = 1.0;
orden = 1;
t = 0.0;

#-

tiempo = @elapsed solfams = Solution_family.(g, x_ini, p_ini, t, Δs, p_fin);

#-

resultado = plot(solfams, leg = false, title = "Familias de Soluciones \n Tiempo de ejecución = $(tiempo) s", ylabel = "u(λ)", xlabel = "λ")
savefig(resultado,"Prueba.png");

#-

precision = []

for i in 1:length(x_ini)
    precisioni = plot(title = "λ_ini = $(p_ini), x_ini = $(x_ini[i])", ylabel = "|g(u(λ),λ)| / 1×10⁻¹⁶", titlefontsize = 7, xlabel = "λ", guidefontsize = 8)
    plot!(solfams[i][1],abs.(g.(solfams[i][2],solfams[i][1],t)) .* 1.e16 ,label = "")
    push!(precision,precisioni)
end

plot( precision..., layout = (2,2),plot_title="Precisión de las ramas de soluciones",plot_titlefontsize = 10)

savefig("Precision.png")