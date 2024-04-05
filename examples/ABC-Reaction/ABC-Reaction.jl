include("../../src/solution_family.jl");
using Plots, Statistics, LinearAlgebra;

#-

function f!(du,u,p,t)
    du[1] = -u[1] + p[1]*(1-u[1])*exp(u[3])
    du[2] = -u[2] + p[1]*(1-u[1])*exp(u[3]) - p[1]*p[4]*u[2]*exp(u[3])
    du[3] = -u[3] - p[3]*u[3] + p[1]*p[5]*(1-u[1])*exp(u[3]) + p[1]*p[5]*p[2]*p[4]*u[2]*exp(u[3])
end;

#-

α = 1;
σ = 0.04;
B = 8;
β = [1.1,1.3,1.5,1.6,1.7,1.8];
D = 0.0;
x_ini = [[0.0,0.0,0.0] for i in β];
p_ini = [[D,α,i,σ,B] for i in β];
Δs= 0.001;
p_fin = 0.5;
indice = 1;
t = 0.0;

#- 

tiempo = @elapsed solfams = Solution_family.(f!, x_ini, p_ini,t, Δs, p_fin,indice);

#-

plot(title = "Reacción A -> B -> C \n TIempo de ejecución = $(tiempo) s", ylabel = "||u(λ)||", xlabel = "λ")

for i in 1:length(solfams)
    plot!(solfams[i][1],norm.(solfams[i][2]), label ="β = $(β[i])")
end

savefig("Prueba.png")

#-

for i in 1:length(x_ini)
    normas = []
    for j in 1:length(solfams[i][1])
        dx = zeros(length(x_ini))
        f!(dx,solfams[i][2][j],[k == indice ? solfams[i][1][j] : p_ini[i][k] for k in 1:length(p_ini[i])],t)
        push!(normas,norm(dx))
    end
    precision = plot(title = "Precisión de la rama", ylabel = "||F(u(λ),λ)||", xlabel = "λ")
    plot!(solfams[i][1],normas ,label = "λ_ini = $(p_ini[i]), x_ini = $(x_ini[i])")
    savefig("PrecisiónRama$(i).png")
end

