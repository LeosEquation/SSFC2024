include("../../src/NumDiffEq.jl");

#-

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
    β = 1.2:0.01:1.42
    D = 0.0
    x_ini = [[0.0,0.0,0.0] for i in β]
    p_ini = [[D,α,i,σ,B] for i in β]
    Δs= 0.001
    p_fin = 0.5
    indice = 1
    t = 0.0
end;

#-

Equilibrium(f!, x_ini[4], p_ini[4],t, Δs, p_fin,indice)

tiempo = @elapsed begin
    ramas_de_equilibrio = Equilibrium.(f!, x_ini, p_ini,t, Δs, p_fin,indice);
    LBP = Limit_Points.(f!,ramas_de_equilibrio,t,p_ini,indice,p_fin);
    HBP = Hopf_Points.(f!,ramas_de_equilibrio,t,p_ini,indice,p_fin);
    estabilidad = Stability_intervals.(f!,p_ini,t,indice,ramas_de_equilibrio);
end;

#-

begin
    LBP = [[(LBP[i][j][1],norm(LBP[i][j][2])) for j in 1:length(LBP[i])] for i in 1:length(LBP)]
    HBP = [[(HBP[i][j][1],norm(HBP[i][j][2])) for j in 1:length(HBP[i])] for i in 1:length(HBP)]
    estable = [[(estabilidad[i][1][j][1],norm(estabilidad[i][1][j][2])) for j in 1:length(estabilidad[i][1])] for i in 1:length(estabilidad)]
    inestable = [[(estabilidad[i][2][j][1],norm(estabilidad[i][2][j][2])) for j in 1:length(estabilidad[i][2])] for i in 1:length(estabilidad)]
end;

#-

begin
    plot(title = "Reacción A -> B -> C \n TIempo de ejecución = $(tiempo) s", ylabel = "||u(D)||", xlabel = "D")
    plot!(estable[1], label = "Estable", linestyle = :solid, color = "blue")
    plot!(inestable[1], label = "Inestable", linestyle = :dash, color = "blue")
    scatter!(HBP[1], label = "Bifurcación Hopf", color = "red")
    scatter!(LBP[1], label = "Bifurcación LP", color = "white")
    for i in 1:length(ramas_de_equilibrio)
        plot!(estable[i], label = "", linestyle = :solid, color = "blue")
        plot!(inestable[i], label = "", linestyle = :dash, color = "blue")
        scatter!(HBP[i], label = "", color = "red")
        scatter!(LBP[i], label = "", color = "white")
    end
    ylims!(1,7)
    xlims!(0.12,0.22)
    savefig("ABC Reaction.png")
end

#-

begin
    #=
    for i in 1:length(x_ini)
        normas = []
        for j in 1:length(ramas_de_equilibrio[i][1])
            dx = zeros(length(x_ini))
            f!(dx,ramas_de_equilibrio[i][2][j],[k == indice ? ramas_de_equilibrio[i][1][j] : p_ini[i][k] for k in 1:length(p_ini[i])],t)
            push!(normas,norm(dx))
        end
        precision = plot(title = "Precisión de la rama", ylabel = "||F(u(λ),λ)||", xlabel = "λ")
        plot!(ramas_de_equilibrio[i][1],normas ,label = "λ_ini = $(p_ini[i]), x_ini = $(x_ini[i])")
        savefig("PrecisiónRama$(i).png")
    end
    =#
end

#-