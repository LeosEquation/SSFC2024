include("../../src/solution_family.jl");
using Plots, Statistics, LinearAlgebra;

#-

function f!(du,u,λ,t)
    du[1] = 3*u[1]*(1 - u[1]) - u[1]*u[2] - λ*(1-exp(-5*u[1]))
    du[2] = -u[2] + 3*u[1]*u[2]
end;

#-

x_ini = [[0.0,0.0],[1.0,0.0],[1/3,2.0]];
p_ini = 0.0;
Δs= 0.001;
p_fin = 0.9;
t = 0.0;

#- 

tiempo = @elapsed solfams = Solution_family.(f!, x_ini, p_ini, t, Δs, p_fin; N = 10000);

#-

plot(title = "Peces \n Tiempo de ejecución = $(tiempo) s", ylabel = "Densidad poblacional de peces", xlabel = "Crecimiento per cápita de la poblacion de peces")
ylims!(-0.1,0.9)
xlims!(0.0,0.9)
for i in 1:length(solfams)
    plot!(solfams[i][1],[j[1] for j in solfams[i][2]], label = "x_ini = $(x_ini[i])")
end
savefig("Prueba.png")

#-

plot(title = "Predator-Prey \n Tiempo de ejecución = $(tiempo) s", ylabel = "Crecimiento per cápita", xlabel = "Densidad poblacional de tiburones", zlabel = "Densidad poblacional de peces", camera = (60, 30))
for i in 1:3
    x = zeros(length(solfams[i][1]))
    y = zeros(length(solfams[i][1]))
    z = zeros(length(solfams[i][1]))

    x = solfams[i][1]

    for j in 1:length(solfams[i][1])
        y[j] = solfams[i][2][j][2]
        z[j] = solfams[i][2][j][1]
    end

    plot!(y,x,z,xflip = true)
end

ylims!(-0.1,1.0)
xlims!(-0.5,2.0)
zlims!(-0.25,1.0)
savefig("Prueba3D.png")

#-

for i in 1:length(x_ini)
    normas = []
    for j in 1:length(solfams[i][1])
        dx = zeros(length(x_ini))
        f!(dx,solfams[i][2][j],solfams[i][1][j],t)
        push!(normas,norm(dx))
    end
    if i == 2
        precision = plot(title = "Precisión de la rama", ylabel = "||F(u(λ),λ)||", xlabel = "λ")
        plot!(solfams[i][1],normas ,label = "λ_ini = $(p_ini), x_ini = $(x_ini[i])")
    else
        precision = plot(title = "Precisión de la rama", ylabel = "||F(u(λ),λ)|| / 1×10⁻¹⁶", xlabel = "λ")
        plot!(solfams[i][1],normas .* 1.e16 ,label = "λ_ini = $(p_ini), x_ini = $(x_ini[i])")
    end
    savefig("PrecisiónRama$(i).png")
end



function f(u,λ,t)
    du1 = 3*u[1]*(1 - u[1]) - u[1]*u[2] - λ*(1-exp(-5*u[1]))
    du2 = -u[2] + 3*u[1]*u[2]
    return [du1,du2]
end;

x = [0.9993710392462667, 1.e-16]    
p = 0.001898527959297181

H = Hessian(f!,x,p,0.0)

S = set_variables("s", numvars = 3, order = 2)

h = TaylorSeries.hessian(sqrt(f(x+S[1:end-1],p+S[end],0.0) ⋅ f(x+S[1:end-1],p+S[end],0.0)))

inv(H)

sqrt(f(x+S[1:end-1],p+S[end],0.0) ⋅ f(x+S[1:end-1],p+S[end],0.0))

4.930380657631324e-32

[2,2] .> 0.0