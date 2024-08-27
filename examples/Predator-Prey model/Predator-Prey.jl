include("../../src/NumDiffEq.jl");

#-

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
    ramas_de_equilibrio = Equilibrium.(f!, x_ini, p_ini, t, Δs, p_fin)
    LBP = Limit_Points.(f!,ramas_de_equilibrio,t,p_fin)
    HBP = Hopf_Points.(f!,ramas_de_equilibrio,t,p_fin)
    estabilidad = Stability_intervals.(f!,t,ramas_de_equilibrio)
    rama_periodica = Periodic_Orbits(f!, HBP[3][1][2], HBP[3][1][1], t, 0.01, p_ini, p_fin; N = 250)
end;

#-

#begin
    plot(title = "Predator-Prey Model \n Tiempo de ejecución = $(tiempo) s", ylabel = "Shark", xlabel = "Fish")
    scatter!(Tuple(HBP[3][1][2]), label = "Bifurcación Hopf", color = "red")
    time = 0.0:rama_periodica[3][2]/100:rama_periodica[3][2]
    sol = taylorinteg(f!, rama_periodica[2][2], time, 20, 1.e-20, rama_periodica[1][2])
    plot!(sol[:,1], sol[:,2], label = "Órbita periódica", color = "blue")

    for i in 3:2:length(rama_periodica[1])
        time = 0.0:rama_periodica[3][i]/100:rama_periodica[3][i]
        sol = taylorinteg(f!, rama_periodica[2][i], time, 20, 1.e-20, rama_periodica[1][i])
        plot!(sol[:,1], sol[:,2], label = "", color = "blue")
    end

    #savefig("Predator-Prey_Periodic_Orbits.png")

    plot!()
#end;

#-

begin
    LBP_2D = [[(LBP[i][j][1],LBP[i][j][2][1]) for j in 1:length(LBP[i])] for i in 1:length(LBP)]
    HBP_2D = [[(HBP[i][j][1],HBP[i][j][2][1]) for j in 1:length(HBP[i])] for i in 1:length(HBP)]
    estable = [[(estabilidad[i][1][j][1],estabilidad[i][1][j][2][1]) for j in 1:length(estabilidad[i][1])] for i in 1:length(estabilidad)]
    inestable = [[(estabilidad[i][2][j][1],estabilidad[i][2][j][2][1]) for j in 1:length(estabilidad[i][2])] for i in 1:length(estabilidad)]
end;



#-

begin
    plot(title = "Predator-Prey Model \n Tiempo de ejecución = $(tiempo) s", ylabel = "Fish", xlabel = "Quota")
    plot!(estable[1], label = "Estable", linestyle = :solid, color = "blue")
    plot!(inestable[1], label = "Inestable", linestyle = :dash, color = "blue")
    scatter!(LBP_2D[1], label = "Bifurcación LP", color = "white")
    scatter!(HBP_2D[1], label = "Bifurcación Hopf", color = "red")
    ylims!(-0.1,1.0)
    xlims!(0.0,0.9)

    for i in 1:length(ramas_de_equilibrio)
        plot!(estable[i], label = "", linestyle = :solid, color = "blue")
        plot!(inestable[i], label = "", linestyle = :dash, color = "blue")
        scatter!(LBP_2D[i], label = "", color = "white")
        scatter!(HBP_2D[i], label = "", color = "red")
    end

    time = [0.0:i/100:i for i in rama_periodica[3]]

    sol = taylorinteg.(f!, rama_periodica[2], time, 20, 1.e-20, rama_periodica[1])

    maxfish = maximum.([i[:,1] for i in sol])
    minfish = minimum.([i[:,1] for i in sol])

    plot!(rama_periodica[1], maxfish, label = "Max Fish")
    plot!(rama_periodica[1], minfish, label = "Min Fish")

    savefig("Predator-Prey_Periodic fish.png")
end;

#-

begin
    LBP_3D = [[(LBP[i][j][2][2],LBP[i][j][1],LBP[i][j][2][1]) for j in 1:length(LBP[i])] for i in 1:length(LBP)]
    HBP_3D = [[(HBP[i][j][2][2],HBP[i][j][1],HBP[i][j][2][1]) for j in 1:length(HBP[i])] for i in 1:length(HBP)]

    plot(title = "Predator-Prey Model \n Tiempo de ejecución = $(tiempo) s", ylabel = "Quota", xlabel = "Sharks", zlabel = "Fish", camera = (60, 25))

    for i in 1:3
        x = zeros(length(ramas_de_equilibrio[i][1]))

        y_estable = zeros(length(ramas_de_equilibrio[i][1]))
        z_estable = zeros(length(ramas_de_equilibrio[i][1]))

        y_inestable = zeros(length(ramas_de_equilibrio[i][1]))
        z_inestable = zeros(length(ramas_de_equilibrio[i][1]))

        x = ramas_de_equilibrio[i][1]

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
            scatter!(LBP_3D[i], label = "Bifurcación LP", color = "white")
            scatter!(HBP_3D[i], label = "Bifurcación Hopf", color = "red")
        else
            plot!(y_estable,x,z_estable,xflip = true, color = "blue",linestyle = :solid, label = "")
            plot!(y_inestable,x,z_inestable,xflip = true, color = "blue",linestyle = :dash, label = "")
            scatter!(LBP_3D[i], label = "", color = "white")
            scatter!(HBP_3D[i], label = "", color = "red")
        end
    end


    ylims!(-0.1,1.0)
    xlims!(-0.5,2.0)
    zlims!(-0.25,1.0)
    savefig("Predator-Prey Model.png")
end

#-

begin
    #=
    for i in 1:length(x_ini)
        normas = []
        for j in 1:length(ramas_de_equilibrio[i][1])
            dx = zeros(length(x_ini))
            f!(dx,ramas_de_equilibrio[i][2][j],ramas_de_equilibrio[i][1][j],t)
            push!(normas,norm(dx))
        end
        if i == 2
            precision = plot(title = "Precisión de la rama", ylabel = "||F(u(λ),λ)||", xlabel = "λ")
            plot!(ramas_de_equilibrio[i][1],normas ,label = "λ_ini = $(p_ini), x_ini = $(x_ini[i])")
        else
            precision = plot(title = "Precisión de la rama", ylabel = "||F(u(λ),λ)|| / 1×10⁻¹⁶", xlabel = "λ")
            plot!(ramas_de_equilibrio[i][1],normas .* 1.e16 ,label = "λ_ini = $(p_ini), x_ini = $(x_ini[i])")
        end
        savefig("PrecisiónRama$(i).png")
    end
    =#
end

#-

#Tiempo de ejecución de 10 órbitas con Taylor1: 41 seg
#Tiempo de ejecución de 10 órbitas con TaylorN: 27 seg
#Tiempo de ejecución de 100 órbitas con TaylorN: 259 seg
#Tiempo de ejecución de 1000 órbitas con TaylorN: 5995 seg