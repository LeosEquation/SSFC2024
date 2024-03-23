include("src/solution_family.jl")
using Plots, Statistics

#-

g(u,λ) = (u^2-1)*(u^2-4) + λ*u^2*exp(u/10)

#-

x_ini = 2.0
p_ini = 0.0
Δs= 0.001
p_fin = 2.0
orden = 1;

#-

@time solfams = Solution_family.(g, x_ini, p_ini, Δs, p_fin);

#-

plot(solfams, leg = false, title = "Familias de Soluciones", ylabel = "u(λ)", xlabel = "λ", c = "blue")

#-

include("src/solution_family.jl")
using Plots, Statistics

#-

g(u,λ) = (u^2-1)*(u^2-4) + λ[1]*u^2*exp(u/10) + sin(λ[2]*u)

#-

x_ini = 2.0
p_ini = [0.0,0.0]
Δs= 0.001
p_fin = 2.0
orden = 1;

#-

@time solfams = Solution_family(g, x_ini, p_ini, Δs, p_fin,1);

#-

plot(solfams, leg = false, title = "Familias de Soluciones", ylabel = "u(λ)", xlabel = "λ", c = "blue")
