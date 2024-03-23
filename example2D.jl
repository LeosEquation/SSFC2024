include("src/solution_family.jl");
using Plots, Statistics, LinearAlgebra;

#-

function f!(du,u,λ)
    du[1] = 3*u[1]*(1 - u[1]) - u[1]*u[2] - λ*(1-exp(-5*u[1]))
    du[2] = -u[2] + 3*u[1]*u[2]
end;

#-

x_ini = [[0.0,0.0],[1.0,0.0],[1/3,2.0]];
p_ini = 0.0;
Δs= 0.001;
p_fin = 0.9;

#- 

@time solfams = Solution_family.(f!, x_ini, p_ini, Δs, p_fin);

#-

plot(title = "Modelo Depredador-Presa", ylabel = "Densidad poblacional de peces", xlabel = "Crecimiento per cápita de la poblacion de peces")
xlims!(0.0,0.9)
ylims!(-0.1,1.2)
for i in 1:length(solfams)
    plot!(solfams[i][1],[j[1] for j in solfams[i][2]], label = "x_ini = $(x_ini[i])")
end
plot!()
#-

