include("../../src/NumDiffEq.jl");

#-

using Plots;

#-

function f!(du,u,λ,t)
    du[1] = 3*u[1]*(1 - u[1]) - u[1]*u[2] - λ*(1-exp(-5*u[1]))
    du[2] = -u[2] + 3*u[1]*u[2]
end;

function df!(J,u,λ,t)
    J[1,1] = 3 - 6*u[1] - u[2] - 5*λ*exp(-5*u[1])
    J[1,2] = -u[1]
    J[2,1] = 3*u[2]
    J[2,2] = -1.0 + 3*u[1]
end;

#-

begin
    x_ini = [[0.0,0.0],[1.0,0.0],[1/3,2.0]]
    p_ini = 0.0
    t = 0.0
    Ne = 5000
    Np = 250
end;

#-

tiempo = @elapsed begin

# Primero calculamos las 3 ramas de equilibrio con sus bifurcaciones de Hopf y LP a un paso de 0.001

ramas_de_equilibrio = Equilibrium.(f!, x_ini, p_ini, t, 0.001, Ne, bifurcations = true);

# Sabemos que en la rama 3 hay una bifurcación de Hopf, entonces calculamos la rama periódica aquí a un paso de 0.01

rama_periodica = Periodic_Orbits(f!, ramas_de_equilibrio[3][8][1, :], ramas_de_equilibrio[3][7][1], t, 0.1, 0.0001, Np; newtontol = 1.e-12, integmaxsteps = 5000)

end;


#-

begin

    plot(title = "Predator-Prey Model \n Tiempo : $(tiempo) s", ylabel = "Fish", xlabel = "Quota")

    p, x, estable, inestable, pf, xf, pb, xb = ramas_de_equilibrio[1]

    plot!([estable[i] ? p[i] : NaN for i in 1:Ne], [estable[i] ? x[i, 1] : NaN for i in 1:Ne], label = "Stable", linestyle = :solid, color = "blue")
    plot!([!estable[i] ? p[i] : NaN for i in 1:Ne], [!estable[i] ? x[i, 1] : NaN for i in 1:Ne], label = "Unstable", linestyle = :dash, color = "blue")

    scatter!(pf, xf[:, 1], label = "Fold bifurcation", color = "white")
    scatter!(pb, xb[:, 1], label = "Hopf bifurcation", color = "red")

    for i in 2:3

        p, x, estable, inestable, pf, xf, pb, xb = ramas_de_equilibrio[i]

        plot!([estable[i] ? p[i] : NaN for i in 1:Ne], [estable[i] ? x[i, 1] : NaN for i in 1:Ne], label = "", linestyle = :solid, color = "blue")
        plot!([!estable[i] ? p[i] : NaN for i in 1:Ne], [!estable[i] ? x[i, 1] : NaN for i in 1:Ne], label = "", linestyle = :dash, color = "blue")

        scatter!(pf, xf[:, 1], label = "", color = "white")
        scatter!(pb, xb[:, 1], label = "", color = "red")

    end

    P, X, T, estable_p = rama_periodica

    times = [0.0:T[i]/1000:T[i] for i in 2:Np]
    sol = taylorinteg.(f!, [X[i,:] for i in 2:Np], times, 20, 1.e-20, P[2:Np])

    plot!(P, [X[1, 1]; [maximum(sol[i][:, 1]) for i in 1:Np-1]], label = "Max órbita periódica")
    plot!(P, [X[1, 1]; [minimum(sol[i][:, 1]) for i in 1:Np-1]], label = "Min órbita periódica")

    ylims!(-0.25, 1.25)
    xlims!(0.0, 0.9)

    savefig("Predator-Prey-fish2.png")

end;

#-

begin

    plot(title = "Predator-Prey Model \n Tiempo : $(tiempo) s", ylabel = "Shark", xlabel = "Fish")

    scatter!((X[1, 1], X[1, 2]), label = "Hopf bifurcation", color = "red")

    plot!(sol[1][:, 1], sol[1][:, 2], label = "Órbitas periódicas", color = "blue")

    for i in 2:Np-1

        plot!(sol[i][:, 1], sol[i][:, 2], label = "", color = "blue")

    end

    savefig("Predator-Prey-periodic2.png")

end;

#-

begin

    plot(title = "Predator-Prey Model \n Tiempo : $(tiempo) s", size=(600, 600), camera = (-30, 15), zlabel = "Fish", ylabel = "Shark", xlabel = "Quota")

    p, x, estable, inestable, pf, xf, pb, xb = ramas_de_equilibrio[1]

    plot!([estable[i] ? p[i] : NaN for i in 1:Ne], [estable[i] ? x[i, 2] : NaN for i in 1:Ne], [estable[i] ? x[i, 1] : NaN for i in 1:Ne], label = "Stable", linestyle = :solid, color = "blue")
    plot!([!estable[i] ? p[i] : NaN for i in 1:Ne], [!estable[i] ? x[i, 2] : NaN for i in 1:Ne], [!estable[i] ? x[i, 1] : NaN for i in 1:Ne], label = "Unstable", linestyle = :dash, color = "blue")

    scatter!(pf, xf[:, 2], xf[:, 1], label = "Fold bifurcation", color = "white")
    scatter!(pb, xb[:, 2], xb[:, 1], label = "Hopf bifurcation", color = "red")

    for i in 2:3

        p, x, estable, inestable, pf, xf, pb, xb = ramas_de_equilibrio[i]

        plot!([estable[i] ? p[i] : NaN for i in 1:Ne], [estable[i] ? x[i, 2] : NaN for i in 1:Ne], [estable[i] ? x[i, 1] : NaN for i in 1:Ne], label = "", linestyle = :solid, color = "blue")
        plot!([!estable[i] ? p[i] : NaN for i in 1:Ne], [!estable[i] ? x[i, 2] : NaN for i in 1:Ne], [!estable[i] ? x[i, 1] : NaN for i in 1:Ne], label = "", linestyle = :dash, color = "blue")

        scatter!(pf, xf[:, 2], xf[:, 1], label = "", color = "white")
        scatter!(pb, xb[:, 2], xb[:, 1], label = "", color = "red")

    end

    P, X, T = rama_periodica

    times = 0.0:T[2]/1000:T[2]
    sol = taylorinteg(f!, X[2,:], times, 20, 1.e-20, P[2])
    plot!([P[2] for i in 1:1001], sol[:, 2], sol[:, 1], label = "Órbita periodica", color = "red")

    for i in 20:20:Np
        times = 0.0:T[i]/1000:T[i]
        sol = taylorinteg(f!, X[i,:], times, 20, 1.e-20, P[i])
        plot!([P[i] for j in 1:1001], sol[:, 2], sol[:, 1], label = "", color = "red")
    end

    xlims!(-0.25, 1.0)
    ylims!(-0.5 , 2.0)
    zlims!(-0.25, 1.0)

    savefig("Predator-Prey Model2.png")

end