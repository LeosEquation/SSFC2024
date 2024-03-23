include("src/solution_family.jl");
using Plots, Statistics, LinearAlgebra;

#-

function f!(du,u,p; α = 1, β = 1.8, σ = 0.04, B = 8)
    du[1] = -u[1] + p*(1-u[1])*exp(u[3])
    du[2] = -u[2] + p*(1-u[1])*exp(u[3]) - p*σ*u[2]*exp(u[3])
    du[3] = -u[3] - β*u[3] + p*B*(1-u[1])*exp(u[3]) + p*B*α*σ*u[2]*exp(u[3])
end;

#-

x_ini = [0.0,0.0,0.0];
p_ini = 0.0;
Δs= 0.001;
p_fin = 0.5;

#- 

@time solfams = Solution_family(f!, x_ini, p_ini, Δs, p_fin);

#-

plot(solfams[1],norm.(solfams[2]), label = "", title = "Familias de Soluciones", ylabel = "||u(λ)||", xlabel = "λ")

#-

function f!(du,u,p)
    du[1] = -u[1] + p[1]*(1-u[1])*exp(u[3])
    du[2] = -u[2] + p[1]*(1-u[1])*exp(u[3]) - p[1]*p[4]*u[2]*exp(u[3])
    du[3] = -u[3] - p[3]*u[3] + p[1]*p[5]*(1-u[1])*exp(u[3]) + p[1]*p[5]*p[2]*p[4]*u[2]*exp(u[3])
end;

#-

α = 1
σ = 0.04
B = 8
β = [1.1,1.3,1.5,1.6,1.7,1.8]
D = 0.0
x_ini = [[0.0,0.0,0.0] for i in β]
p_ini = [[D,α,i,σ,B] for i in β]
Δs= 0.001
p_fin = 0.5
indice = 1;

#- 

@time solfams = Solution_family.(f!, x_ini, p_ini, Δs, p_fin,indice);

#-

plot(title = "Familias de Soluciones", ylabel = "||u(λ)||", xlabel = "λ")

for i in 1:length(solfams)
    plot!(solfams[i][1],norm.(solfams[i][2]), label ="β = $(β[i])")
end

plot!()

#-