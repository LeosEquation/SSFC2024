using TaylorSeries, LinearAlgebra, Plots

include("src/solution_family.jl")

function f!(du,u,λ,t)
    du[1] = 3*u[1]*(1 - u[1]) - u[1]*u[2] - λ*(1-exp(-5*u[1]))
    du[2] = -u[2] + 3*u[1]*u[2]
end;

function f(u,λ,t)
    du1 = 3*u[1]*(1 - u[1]) - u[1]*u[2] - λ*(1-exp(-5*u[1]))
    du2 = -u[2] + 3*u[1]*u[2]

    return [du1,du2]
end;

x_ini = [1.0,0.0];
p_ini = 0.0;
Δs= 0.001;
p_fin = 0.9;
t = 0.0;

#-

solfams = Solution_family(f!, x_ini, p_ini, t, Δs, p_fin);

#-

J = Jacobian.(f!,solfams[2],solfams[1],t)
J = [j[:,1:end-1] for j in J]
realparts = real.(eigvals.(J))

realparts_p = [i[1] for i in realparts]
realparts_n = [i[2] for i in realparts]

plot(solfams[1],realparts_p)
plot!(solfams[1],realparts_n)

imagparts = imag.(eigvals.(J))

imagparts_p = [i[1] for i in imagparts]
imagparts_n = [i[2] for i in imagparts]

plot!(solfams[1],imagparts_p)
plot!(solfams[1],imagparts_n)

ylims!(-5,5)


findmin(abs.(realparts_p))
findmin(abs.(realparts_n))

realparts_n[1557]

solfams[1][1557]

plot(title = "Peces \n Tiempo de ejecución = $(tiempo) s", ylabel = "Densidad poblacional de peces", xlabel = "Crecimiento per cápita de la poblacion de peces")
ylims!(-0.1,0.9)
xlims!(0.0,0.9)
for i in 1:length(solfams)
    plot!(solfams[i][1],[j[1] for j in solfams[i][2]], label = "x_ini = $(x_ini[i])")
end
vline!([solfams[3][1][1768]])

minimos = []
indices = []

for i in 3:(length(solfams[1])-2)
    if abs.(realparts_n)[i-1] > abs.(realparts_n)[i] < abs.(realparts_n)[i+1]
        push!(minimos,abs.(realparts_n)[i])
        push!(indices,i)
    end
end

minimos

indices

plot(solfams[1],norm.(solfams[2]))
vline!(solfams[1][indices])

solfams[1][indices]

dx = x_ini
normas = []
for i in 1:length(solfams[1])
    f!(dx,solfams[2][i],[j == indice ? solfams[1][i] : p_ini[j] for j in 1:length(p_ini)],t)
    push!(normas,norm(dx))
end

plot(solfams[1],normas)

findmax(norm.(f.(solfams[2],solfams[1],t)))
norm.(f.(solfams[2],solfams[1],t))[7022]

plot(-5.2967:0.000001:-5.2956,-0.005:0.000001:0.0,g,st=:surface,aspect_ratio=:equal)
scatter!((solfams[2][7022][1],solfams[2][7022][2],0.0))
zlims!(0.0,7.685230256981868e-5)

norm.(f.(solfams[1],solfams[2],t))

g(x,y) = norm(f([x,y],solfams[1][7022],t))^2

g(solfams[2][7022]...)

solfams[2][7022]

plot(solfams[1],norm.(f.(solfams[2],solfams[1],t)))

[solfams[2][7022]; 0.0]