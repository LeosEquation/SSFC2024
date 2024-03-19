include("src/solution_family.jl")
include("src/bifurcation.jl")
using Plots

#-

g(u,λ) = (u^2-1)*(u^2-4) + λ*u^2*exp(u/10)

#-

x_ini = [1.0,-1.0]
p_ini = 0.0
Δs= 0.001
p_fin = 2.0
orden = 10;

#-

@time solfams = Solution_family.(g, x_ini, p_ini, Δs, p_fin,orden);

#-

plot(solfams, leg = false, title = "Familias de Soluciones", ylabel = "u(λ)", xlabel = "λ", c = "blue")

#-

Bifurcation_point(g,x_ini[1],p_ini,orden)

#_

function f!(du,u,p)
    du[1] = -u[1] + p[1]*(1-u[1])*exp(u[3])
    du[2] = -u[2] + p[1]*(1-u[1])*exp(u[3]) - p[1]*p[4]*u[2]*exp(u[3])
    du[3] = -u[3] - p[3]*u[3] + p[1]*p[5]*(1-u[1])*exp(u[3]) + p[1]*p[5]*p[2]*p[4]*u[2]*exp(u[3])
end

#-

α = 1
σ = 0.04
B = 8
D = 0.0
x_ini = [0.0,0.0,0.0]
Δs= 0.001
p_fin = 0.5
orden = 10
indice = 1;

#- 

plot(title = "Familias de Soluciones", ylabel = "||u(λ)||", xlabel = "λ")

@time for β in [1.1,1.3,1.5,1.6,1.7,1.8]
    p_ini = [D,α,β,σ,B]
    solfams = Solution_family(f!, x_ini, p_ini, Δs, p_fin,orden,indice);
    plot!(solfams[1],norm.(solfams[2]), label = "β = $(β)")
end

plot!()

#-

using LinearAlgebra
p_ini = [D,α,1.1,σ,B]
solfams = Solution_family(f!, x_ini, p_ini, Δs, p_fin,orden,indice)
normas = Float64[]
for i in 1:length(solfams[1])
    dx = [0.0,0.0,0.0]
    f!(dx,solfams[2][i],[solfams[1][i],α,1.1,σ,B])
    push!(normas,norm(dx))
end

#-

plot(normas)