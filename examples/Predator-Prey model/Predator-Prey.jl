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
end;

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

    savefig("Predator-Prey fish.png")
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

prueba = Equilibrium(f!,x_ini[3],p_ini,t,Δs,p_fin)
pb, xb = Hopf_Points(f!,prueba,t,p_fin)[1]

#-

function fp!(du,u,λ,t)
    du[1] = - 1 + exp(-5*u[1])
    du[2] = 0.0
end;

#-

function G(f!,x,p,T,Δs,x0,p0,T0; orden = 20, tol = 1.e-20)
    Time, X = integracion_taylor(f!,x,0.0,T,20,1.e-20,p; Nt = 100)
    Time0, X0 = integracion_taylor(f!,x0,0.0,T,20,1.e-20,p0; Nt = 100)

    G1 = X[end] - x

    dX = [zeros(length(i)) for i in X]
    dX0 = [zeros(length(i)) for i in X0]

    for i in 1:length(Time)
        dx = zeros(length(x))
        f!(dx,X[i],p,Time[i])
        dX[i] = copy(dx)
    end

    for i in 1:length(Time0)
        dx = zeros(length(x))
        f!(dx,X0[i],p0,Time0[i])
        dX0[i] = copy(dx)
    end

    G2 = 0.0

    for i in 2:length(Time)

        G2 += 0.5*(Time[i]*Time[i-1])*((X[i-1] ⋅ dX0[i-1]) + (X[i] ⋅ dX0[i]))

    end

    #G3 = (x-x0) ⋅ x_s + (T - T0)*T_s + (p-p0)*p_s - Δs

    return [G1;G2]

end

#-

G(f!,xb,pb + 0.001 ,T + 0.001,0.001,xb,pb,T)

#-

n = length(xb)
Jx = zeros(n,n)
Jp = zeros(n)
F0 = zeros(n)
s = Taylor1(2)
r = Taylor1([0.0,0.0],2)

#-

x = xb
p = pb
t = 0.0

#-

for i in 1:n
    for j in 1:n
        dx = [s for i in 1:n]
        f!(dx,x+[j == k ? s : r for k in 1:n],p + r, t + r)
        Jx[i,j] = derivative(dx[i])(0.0)
    end
end

#-

for i in 1:n
    dx = [s for i in 1:n]
    f!(dx,x .+ r,p + s,t + r)
    Jp[i] = derivative(dx[i])(0.0)
end

#-

T = 2*π / abs(imag(eigvals(Jx)[1]))

#-

h = (1.e-16)^(1/3) * norm(pb)

#-

x = [1.0,1.0]
p = 1.0
h = 1.e-12

#-

Time1, X = integracion_taylor(f!,x,0.0,T,20,1.e-20,p; Nt = 100)
Timeb, Xb = integracion_taylor(f!,x,0.0,T,20,1.e-20,p - h; Nt = 100)
Timef, Xf = integracion_taylor(f!,x,0.0,T,20,1.e-20,p + h; Nt = 100)
Time2, X2 = integracion_taylor(fp!,x,0.0,T,20,1.e-20,p; Nt = 100)
Time3, X3 = integracion_taylor(fp1!,x,0.0,T,20,1.e-20,p; Nt = 100)

#-

Xp = (Xf .- Xb) ./ (2*h)

#-

plot(Time, [[i[1] for i in X],[i[1] for i in Xb],[i[1] for i in Xf]])

#-

plot(Time, [[i[2] for i in X],[i[2] for i in Xb],[i[2] for i in Xf]])

#-

plot(Time,[[i[1]-x[1] for i in X2],[i[1] for i in Xp]])

#-

function g!(dx,x,p,t)
    f!(dx,x,p,t)
    dx[1] = T*dx[1]
    dx[2] = T*dx[2]
end

#-