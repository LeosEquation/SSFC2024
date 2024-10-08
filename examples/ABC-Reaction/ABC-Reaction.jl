include("../../src/NumDiffEq.jl");


#-

using Plots;

#-

function f!(dx,x,p,t)
    dx[1] = -x[1] + p[1]*(1-x[1])*exp(x[3])
    dx[2] = -x[2] + p[1]*(1-x[1])*exp(x[3]) - p[1]*p[4]*x[2]*exp(x[3])
    dx[3] = -x[3] - p[3]*x[3] + p[1]*p[5]*(1-x[1])*exp(x[3]) + p[1]*p[5]*p[2]*p[4]*x[2]*exp(x[3])
    return dx
end;

#-

begin
    α = 1
    σ = 0.04
    B = 8
    β = 1.55
    D = 0.0
    x_ini = [0.0,0.0,0.0] # [[0.0,0.0,0.0] for i in β]
    p_ini = [D,α,β,σ,B] #[[D,α,i,σ,B] for i in β]
    Δs= 0.001
    p_fin = 0.5
    indice = 1
    N = 6000
    t = 0.0
    Np = 50
    Δsp = 0.05
    Δsp_ini = 1.e-3
end;

#-

@elapsed begin

@time p, x, estable, inestable, pf, xf, pb, xb = Equilibrium(f!, x_ini, p_ini, t, Δs, indice, N; bifurcations = true);

@time P1, X1, T1, stable1 = Periodic_Orbits(f!, xb[1,:], pb[1,:], t, 0.01, 0.00001, indice, 500; newtontol = 1.e-12, integmaxsteps = 2000);

@time P2, X2, T2, stable2 = Periodic_Orbits(f!, xb[2,:], pb[2,:], t, 0.01, 0.00001, indice, 500; newtontol = 1.e-12, integmaxsteps = 2000);

@time P3, X3, T3, stable3 = Periodic_Orbits(f!, xb[3,:], pb[3,:], t, 0.01, 0.00001, indice, 500; newtontol = 1.e-12, integmaxsteps = 2000);

@time P4, X4, T4, stable4 = Periodic_Orbits(f!, xb[4,:], pb[4,:], t, 0.01, 0.00001, indice, 500; newtontol = 1.e-12, integmaxsteps = 2000);

end;

#-

times1 = [0.0:i/1000:i for i in T1];
sol1 = taylorinteg.(f!,[X1[i,:] for i in 1:length(P1[:, 1])] , times1, 20, 1.e-20, [P1[i,:] for i in 1:length(P1[:, 1])]; parse_eqs = false);

times2 = [0.0:i/1000:i for i in T2];
sol2 = taylorinteg.(f!,[X2[i,:] for i in 1:length(P2[:, 1])] , times2, 20, 1.e-20, [P2[i,:] for i in 1:length(P2[:, 1])]);

times3 = [0.0:i/1000:i for i in T3];
sol3 = taylorinteg.(f!,[X3[i,:] for i in 1:length(P3[:, 1])] , times3, 20, 1.e-20, [P3[i,:] for i in 1:length(P3[:, 1])]);

times4 = [0.0:i/1000:i for i in T4];
sol4 = taylorinteg.(f!,[X4[i,:] for i in 1:length(P4[:, 1])] , times4, 20, 1.e-20, [P4[i,:] for i in 1:length(P4[:, 1])]);

#-

i = 2
j = 3

scatter((xb[:,i], xb[:,j]), label = "Bifurcación Hopf", color = "red")

plot!(sol1[2][:,i], sol1[2][:,j], label = "Órbita periódica", color = "blue")

for ii in 50:50:length(T1)
    plot!(sol1[ii][:,i], sol1[ii][:,j], label = "", color = "blue")
end
plot!(X1[:,i], X1[:,j], color = "red", label = "x(0)")

plot!(sol2[2][:,i], sol2[2][:,j], label = "Órbita periódica", color = "blue")

for ii in 50:50:length(T2)
    plot!(sol2[ii][:,i], sol2[ii][:,j], label = "", color = "blue")
end

plot!(sol3[2][:,i], sol3[2][:,j], label = "Órbita periódica", color = "blue")

for ii in 50:50:length(T3)
    plot!(sol3[ii][:,i], sol3[ii][:,j], label = "", color = "blue")
end

plot!(sol4[2][:,i], sol4[2][:,j], label = "Órbita periódica", color = "blue")


for ii in 50:50:length(T4)
    plot!(sol4[ii][:,i], sol4[ii][:,j], label = "", color = "blue")
end

plot!(X2[:,i], X2[:,j], color = "red", label = "x(0)")
plot!(X3[:,i], X3[:,j], color = "red", label = "x(0)")
plot!(X4[:,i], X4[:,j], color = "red", label = "x(0)")

plot!(title = "ABC Reaction", xlabel = "x$(i)", ylabel = "x$(j)",size = (1000,1000))

savefig("ABC$(i)$(j).png")

#-

plot(title = "ABC Reaction", xlabel = "x1", ylabel = "x2", zlabel = "x3", size = (1000,1000))

plot!(sol1[2][:,1], sol1[2][:,2], sol1[2][:,3], label = "Órbita periódica", color = "blue")
plot!(sol2[2][:,1], sol2[2][:,2], sol2[2][:,3], label = "Órbita periódica", color = "green")
plot!(sol3[2][:,1], sol3[2][:,2], sol3[2][:,3], label = "Órbita periódica", color = "black")
plot!(sol4[2][:,1], sol4[2][:,2], sol1[2][:,3], label = "Órbita periódica", color = "orange")

for ii in 50:50:length(P1[:,1])
    plot!(sol1[ii][:,1], sol1[ii][:,2],sol1[ii][:,3], label = "", color = "blue")
end

for ii in 50:50:length(P2[:,1])
    plot!(sol2[ii][:,1], sol2[ii][:,2],sol2[ii][:,3], label = "", color = "green")
end

for ii in 50:50:length(P3[:,1])
    plot!(sol3[ii][:,1], sol3[ii][:,2],sol3[ii][:,3], label = "", color = "black")
end

for ii in 50:50:length(P4[:,1])
    plot!(sol4[ii][:,1], sol4[ii][:,2],sol4[ii][:,3], label = "", color = "orange")
end

scatter!(xb[:, 1], xb[:, 2], xb[:, 3], label = "Bifurcación Hopf", color = "red")

savefig("ABC3d.png")

#-

plot(title = "ABC Reaction \n  s", xlabel = "D", ylabel = "x3", size = (1200, 1200))

plot!([estable[i] ? p[i, indice] : NaN for i in 1:N], [estable[i] ? x[i, 3] : NaN for i in 1:N], label = "Stable", linestyle = :solid, color = "blue")
plot!([!estable[i] ? p[i, indice] : NaN for i in 1:N], [!estable[i] ? x[i, 3] : NaN for i in 1:N], label = "Unstable", linestyle = :dash, color = "blue")

scatter!(pf[:, indice], xf[:, 3], label = "Fold bifurcation", color = "white")
scatter!(pb[:, indice], xb[:, 3], label = "Hopf bifurcation", color = "red")

plot!(P1[:, indice], [X1[1, 3] ; maximum.([sol1[i][:, 3] for i in 2:length(T1)])], label = "Máximo órbita periódica", color = "red")

plot!(P2[:, indice], [X2[1, 3] ; maximum.([sol2[i][:, 3] for i in 2:length(T2)])], label = "", color = "red")

plot!(P3[:, indice], [X3[1, 3] ; maximum.([sol3[i][:, 3] for i in 2:length(T3)])], label = "", color = "red")

plot!(P4[:, indice], [X4[1, 3] ; maximum.([sol4[i][:, 3] for i in 2:length(T4)])], label = "", color = "red")

plot!(P1[:, indice], [X1[1, 3] ; minimum.([sol1[i][:, 3] for i in 2:length(T1)])], label = "Mínimo órbita periódica", color = "orange")

plot!(P2[:, indice], [X2[1, 3] ; minimum.([sol2[i][:, 3] for i in 2:length(T2)])], label = "", color = "orange")

plot!(P3[:, indice], [X3[1, 3] ; minimum.([sol3[i][:, 3] for i in 2:length(T3)])], label = "", color = "orange")

plot!(P4[:, indice], [X4[1, 3] ; minimum.([sol4[i][:, 3] for i in 2:length(T4)])], label = "", color = "orange")

ylims!(1, 10)
xlims!(0.2, 0.4)

savefig("ABC4.png")

#-

plot(P1[:, 1], [norm(X1[i, :]) for i in 1:length(T1)], ylabel = "||x(0)||", xlabel = "p", leg = false)
scatter!(P1[:, 1], [norm(X1[i, :]) for i in 1:length(T1)], ylabel = "||x(0)||", xlabel = "p", leg = false, markersize = 2.0)
savefig("ABCInitial1.png")

plot(P2[:, 1], [norm(X2[i, :]) for i in 1:length(T2)], ylabel = "||x(0)||", xlabel = "p", leg = false)
scatter!(P2[:, 1], [norm(X2[i, :]) for i in 1:length(T2)], ylabel = "||x(0)||", xlabel = "p", leg = false, markersize = 2.0)
savefig("ABCInitial2.png")

plot(P3[:, 1], [norm(X3[i, :]) for i in 1:length(T3)], ylabel = "||x(0)||", xlabel = "p", leg = false)
scatter!(P3[:, 1], [norm(X3[i, :]) for i in 1:length(T3)], ylabel = "||x(0)||", xlabel = "p", leg = false, markersize = 2.0)
savefig("ABCInitial3.png")

plot(P4[:, 1], [norm(X4[i, :]) for i in 1:length(T4)], ylabel = "||x(0)||", xlabel = "p", leg = false)
scatter!(P4[:, 1], [norm(X4[i, :]) for i in 1:length(T4)], ylabel = "||x(0)||", xlabel = "p", leg = false, markersize = 2.0)
savefig("ABCInitial4.png")