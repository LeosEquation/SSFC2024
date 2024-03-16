include("src/solution_family.jl")
using .SolutionFamily, Plots

g(u,λ) = (u^2-1)*(u^2-4) + λ*u^2*exp(u/10)

x_ini = [1.0,-1.0]
p_ini = 0.0
Δs= 0.001
p_fin = 2.0
orden = 10;

@time solfams = Solution_family.(g, x_ini, p_ini, Δs, p_fin,orden);

plot(solfams, leg = false, title = "Familias de Soluciones", ylabel = "u(λ)", xlabel = "λ", c = "blue")

function f!(du,u,λ)
    du[1] = 3*u[1]*(1-u[1]) - u[1]*u[2] - λ*(1-exp(-5*u[1]))
    du[2] = -u[2] + 3*u[1]*u[2]
end

x_ini = [0.0,0.0]
p_ini = 0.0
Δs= 0.001
p_fin = 1.0
orden = 10;

@time solfams = Solution_family(f!, x_ini, p_ini, Δs, p_fin,orden);

transpose([1,2,3])*[4,5,6]