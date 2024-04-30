include("../../src/equilibrium.jl");
include("../../src/bifurcation.jl");
include("../../src/stability.jl");

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
    β = 1.2:0.01:1.4
    D = 0.0
    x_ini = [[0.0,0.0,0.0] for i in β]
    p_ini = [[D,α,i,σ,B] for i in β]
    Δs= 0.001
    p_fin = 0.22
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
    plot(title = "Reacción A -> B -> C \n TIempo de ejecución = $(tiempo) s", ylabel = "||u(D)||", xlabel = "D")
    plot!(estable[1], label = "Estable", linestyle = :solid, color = "blue")
    plot!(inestable[1], label = "Inestable", linestyle = :dash, color = "blue")
    scatter!(Hb[1], label = "Bifurcación de Hopf", color = "red")
    for i in 1:length(solfams)
        plot!(estable[i], label = "", linestyle = :solid, color = "blue")
        plot!(inestable[i], label = "", linestyle = :dash, color = "blue")
        scatter!(Hb[i], label = "", color = "red")
    end
    plot!()
end

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

#-

function BiAlternate(A)
    n = size(A)[1]
    m = Int64((n-1)*n/2)
    M = [A[1,1] for i in 1:m, j in 1:m]
    indices = [Int64[] for i in 1:m]

    i = 1
    for j in 2:n
        for k in 1:(j-1)
            indices[i] = [k,j]
            i += 1
        end
    end

    Indices = [[i,j] for i in indices, j in indices]

    for i in 1:m
        for j in 1:m

            p = Indices[i,j][1][1]
            q = Indices[i,j][1][2]
            r = Indices[i,j][2][1]
            s = Indices[i,j][2][2]

            if r == q
                M[i,j] = -A[p,s]
            elseif r != p && s == q
                M[i,j] = A[p,r]
            elseif r == p && s == q
                M[i,j] = A[p,p] + A[q,q]
            elseif r == p && s != q
                M[i,j] = A[q,s]
            elseif s == p
                M[i,j] = -A[q,r]
            else
                M[i,j] = A[p,q]*0.0
            end
        end
    end

    return M

end

S = set_variables("x", numvars = 3, order = 10)

A = [(i - j)^2 + S[i] for i in 1:3,j in 1:3]

BA = BiAlternate(A)

Determinant(A) = A[1,1]*A[2,2]*A[3,3] + A[1,2]*A[2,3]*A[1,3] + A[2,1]*A[3,2]*A[3,1] - A[1,3]*A[2,2]*A[3,1] - A[1,2]*A[2,1]*A[3,3] - A[2,3]*A[3,2]*A[1,1]

Determinant(BA)

det(BA(zeros(3)))

prueba = Equilibrium(f!,x_ini[1],p_ini[1],t,Δs,p_fin,indice)

estabilidad = Stability_intervals(f!,p_ini[1],t,indice,prueba)

begin
    estable = [(estabilidad[1][j][1],norm(estabilidad[1][j][2])) for j in 1:length(estabilidad[1])]
    inestable = [(estabilidad[2][j][1],norm(estabilidad[2][j][2])) for j in 1:length(estabilidad[2])]
end;


begin
    plot(title = "Reacción A -> B -> C", ylabel = "||u(D)||", xlabel = "D")
    plot!(estable, label = "Estable", linestyle = :solid, color = "blue")
    plot!(inestable, label = "Inestable", linestyle = :dash, color = "blue")
end

function Bif_vec(f!,x0,p0,t,indice)
    S = set_variables("x", numvars = length(x0)+1, order = 3)
    x = x0 + S[1:end-1]
    p = p0 + [i == indice ? S[end] : 0.0 for i in 1:length(p0)]
    dx = x
    f!(dx,x,p,t)
    J = [derivative(dx[i],j) for i in 1:length(dx), j in 1:length(dx)]

    return [dx;Determinant(BiAlternate(J))]
end

Hb = Hopf_bifurcation(f!,p_ini[1],t,indice,prueba)


F = Bif_vec(f!,Hb[1][2],[i == indice ? Hb[1][1] : p_ini[1][i] for i in 1:length(p_ini[1])],t,indice)

i = 1

inv(TaylorSeries.jacobian(F))

x = Hb[1][2]
p = [i == indice ? Hb[1][1] : p_ini[1][i] for i in 1:length(p_ini[1])]

while i <= 30
    NewtonVec = [x;p[indice]] - inv(TaylorSeries.jacobian(F))*F(zeros(length(x)+1))
    x = NewtonVec[1:end-1]
    p = [i == indice ? NewtonVec[end] : p_ini[1][i] for i in 1:length(p_ini[1])]
    F = Bif_vec(f!,x,p,t,indice)
    println(norm(F(zeros(length(x) + 1))))
    i += 1
end