include("../../src/equilibrium.jl");
include("../../src/bifurcation.jl");
include("../../src/stability.jl");
using Plots;

#-

function f!(du,u,p,t)
    du[1] = -u[1] + p[1]*(1-u[1])*exp(u[3])
    du[2] = -u[2] + p[1]*(1-u[1])*exp(u[3]) - p[1]*p[4]*u[2]*exp(u[3])
    du[3] = -u[3] - p[3]*u[3] + p[1]*p[5]*(1-u[1])*exp(u[3]) + p[1]*p[5]*p[2]*p[4]*u[2]*exp(u[3])
end;

#-
begin
    α = 1
    σ = 0.04
    B = 8
    β = [1.1,1.3,1.5,1.6,1.7,1.8]
    D = 0.0
    x_ini = [[0.0,0.0,0.0] for i in β]
    p_ini = [[D,α,i,σ,B] for i in β]
    Δs= 0.001
    p_fin = 0.5
    indice = 1
    t = 0.0
end;

#- 

tiempo = @elapsed begin
    solfams = Equilibrium.(f!, x_ini, p_ini,t, Δs, p_fin,indice);
    Hb = Hopf_bifurcation.(f!,p_ini,t,indice,solfams)
    estabilidad = Stability_intervals.(f!,p_ini,t,indice,solfams) 
end;

#-

begin
    Hb = [[(Hb[i][j][1],norm(Hb[i][j][2])) for j in 1:length(Hb[i])] for i in 1:length(Hb)]
    estable = [[(estabilidad[i][1][j][1],norm(estabilidad[i][1][j][2])) for j in 1:length(estabilidad[i][1])] for i in 1:length(estabilidad)]
    inestable = [[(estabilidad[i][2][j][1],norm(estabilidad[i][2][j][2])) for j in 1:length(estabilidad[i][2])] for i in 1:length(estabilidad)]
end;

#-

begin
    plot(title = "Reacción A -> B -> C \n TIempo de ejecución = $(tiempo) s", ylabel = "||u(λ)||", xlabel = "λ")
    plot!(estable[1], label = "Estable", linestyle = :solid, color = "blue")
    plot!(inestable[1], label = "Inestable", linestyle = :dash, color = "blue")
    scatter!(Hb[1], label = "Bifurcación de Hopf", color = "red")
    for i in 1:length(solfams)
        plot!(estable[i], label = "", linestyle = :solid, color = "blue")
        plot!(inestable[i], label = "", linestyle = :dash, color = "blue")
        scatter!(Hb[i], label = "", color = "red")
    end

    savefig("ABC Reaction.png")
end;

#-

begin
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
end
