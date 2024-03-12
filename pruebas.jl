include("psarc.jl")
using .PArcLength, Plots

g(u,λ) = (u^2-1)*(u^2-4) + λ*u^2*exp(u/10)

x_ini = [2.0,1.0,-1.0,-2.0]
p_ini = 0.0
Δs= 0.001
p_fin = 1.0
orden = 10

solfams = SolFam.(g, x_ini, p_ini, Δs, p_fin,orden)

plot(solfams, leg = false, title = "Familias de Soluciones", ylabel = "u(λ)", xlabel = "λ", c = "blue")