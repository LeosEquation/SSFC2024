include("../../src/equilibrium.jl");
include("../../src/bifurcation.jl");
include("../../src/stability.jl");
using Plots;

#-

function f!(du,u,λ,t)
    du[1] = 3*u[1]*(1 - u[1]) - u[1]*u[2] - λ*(1-exp(-5*u[1]))
    du[2] = -u[2] + 3*u[1]*u[2]
end;

#-

begin
    x_ini = [[0.0,0.0],[1.0,0.0],[1/3,2.0]]
    p_ini = 0.0
    Δs= 0.001
    p_fin = 0.9
    t = 0.0
end;

#-

tiempo = @elapsed begin
    solfams = Equilibrium.(f!, x_ini, p_ini, t, Δs, p_fin)
    Hb = Hopf_bifurcation.(f!,t,solfams)
    estabilidad = Stability_intervals.(f!,t,solfams)
end;

#-

begin
    Hb2D = [[(Hb[i][j][1],Hb[i][j][2][1]) for j in 1:length(Hb[i])] for i in 1:length(Hb)]
    estable = [[(estabilidad[i][1][j][1],estabilidad[i][1][j][2][1]) for j in 1:length(estabilidad[i][1])] for i in 1:length(estabilidad)]
    inestable = [[(estabilidad[i][2][j][1],estabilidad[i][2][j][2][1]) for j in 1:length(estabilidad[i][2])] for i in 1:length(estabilidad)]
end;

#-

begin
    plot(title = "Predator-Prey Model \n Tiempo de ejecución = $(tiempo) s", ylabel = "Fish", xlabel = "Quota")
    plot!(estable[1], label = "Estable", linestyle = :solid, color = "blue")
    plot!(inestable[1], label = "Inestable", linestyle = :dash, color = "blue")
    scatter!(Hb2D[1], label = "Bifurcación de Hopf", color = "red")
    ylims!(-0.1,1.0)
    xlims!(0.0,0.9)

    for i in 1:length(solfams)
    plot!(estable[i], label = "", linestyle = :solid, color = "blue")
    plot!(inestable[i], label = "", linestyle = :dash, color = "blue")
    scatter!(Hb2D[i], label = "", color = "red")
    end

    savefig("Predator-Prey fish.png")
end;

#-

begin
    Hb3D = [[(Hb[i][j][2][2],Hb[i][j][1],Hb[i][j][2][1]) for j in 1:length(Hb[i])] for i in 1:length(Hb)]

    plot(title = "Predator-Prey Model \n Tiempo de ejecución = $(tiempo) s", ylabel = "Quota", xlabel = "Sharks", zlabel = "Fish", camera = (60, 25))

    for i in 1:3
        x = zeros(length(solfams[i][1]))

        y_estable = zeros(length(solfams[i][1]))
        z_estable = zeros(length(solfams[i][1]))

        y_inestable = zeros(length(solfams[i][1]))
        z_inestable = zeros(length(solfams[i][1]))

        x = solfams[i][1]

        for j in 1:length(estabilidad[i][1])
            y_estable[j] = estabilidad[i][1][j][2][2]
            z_estable[j] = estabilidad[i][1][j][2][1]
        end

        for j in 1:length(estabilidad[i][2])
            y_inestable[j] = estabilidad[i][2][j][2][2]
            z_inestable[j] = estabilidad[i][2][j][2][1]
        end

        if i == 1
            plot!(y_estable,x,z_estable,xflip = true, color = "blue",linestyle = :solid, label = "Estable")
            plot!(y_inestable,x,z_inestable,xflip = true, color = "blue",linestyle = :dash, label = "Inestable")
            scatter!(Hb3D[i], label = "Bifurcación de Hopf", color = "red")
        else
            plot!(y_estable,x,z_estable,xflip = true, color = "blue",linestyle = :solid, label = "")
            plot!(y_inestable,x,z_inestable,xflip = true, color = "blue",linestyle = :dash, label = "")
            scatter!(Hb3D[i], label = "", color = "red")
        end
    end


    ylims!(-0.1,1.0)
    xlims!(-0.5,2.0)
    zlims!(-0.25,1.0)
    savefig("Predator-Prey Model.png")
end

#-

begin
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
end
