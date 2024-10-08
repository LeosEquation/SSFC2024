# # Bifurcaciones
# Este archivo contiene las funciones correspondientes para encontrar
# puntos de bifurcación en un sistema de ecuaciones diferenciales.

#-

function Limit_Points!(f::Function, pf::Vector{Float64}, xf::Vector{Float64}; newtonite = 8, newtontol = 1.e-16)

    m = length(xf)

    w = 1.0

    F = zeros(3)
    J = zeros(3, 3)

    S = set_variables("s", numvars = 2, order = 2)

    r = zero(S[1])

    dfx = 0.0
    dfp = 0.0
    dfxx = 0.0
    dfxp = 0.0

    zerosTN = zeros(2)

    Δ = zeros(3)

    if m != 0

        for i in 1:m

            println(" \n Fold $(i) : \n ")

            dx = f(xf[i] + S[1], pf[i] + S[2], t + r)

            dfx = differentiate(dx, 1)(zerosTN)
            dfp = differentiate(dx, 2)(zerosTN)
            dfxx = differentiate(differentiate(dx, 1),1)(zerosTN)
            dfxp = differentiate(differentiate(dx, 1),2)(zerosTN)

            v = dfx

            LP_Function!(F, f, dfx, xf[i], pf[i], v, w, t)
            LP_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w)

            j = 1

            while j <= newtonite && norm(F) > newtontol

                Δ .= - inv(J)*F

                xf[i] += Δ[1]
                pf[i] += Δ[2]
                v += Δ[3]

                dx = f(xf[i] + S[1], pf[i] + S[2], t + r)

                dfx = differentiate(dx, 1)(zerosTN)
                dfp = differentiate(dx, 2)(zerosTN)
                dfxx = differentiate(differentiate(dx, 1),1)(zerosTN)
                dfxp = differentiate(differentiate(dx, 1),2)(zerosTN)

                LP_Function!(F, f, dfx, xf[i], pf[i], v, w, t)
                LP_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w)

                j += 1

                println("$(j) : $(norm(F))")

            end

        end

    end
    
end

#-

function Limit_Points!(f::Function, pf::Matrix{Float64}, xf::Vector{Float64}, indice::Int64; newtonite = 8, newtontol = 1.e-16)

    m = length(xf)
    mp, np = size(pf)

    w = 1.0

    F = zeros(3)
    J = zeros(3, 3)

    S = set_variables("s", numvars = 2, order = 2)
    r = zero(S[1])

    Sp = [i == indice ? S[end] : r for i in 1:np]

    dfx = 0.0
    dfp = 0.0
    dfxx = 0.0
    dfxp = 0.0

    zerosTN = zeros(2)

    Δ = zeros(3)

    if m != 0

        for i in 1:m

            println(" \n Fold $(i) : \n ")

            dx = f(xf[i] + S[1], pf[i, :] .+ Sp, t + r)

            dfx = differentiate(dx, 1)(zerosTN)
            dfp = differentiate(dx, 2)(zerosTN)
            dfxx = differentiate(differentiate(dx, 1),1)(zerosTN)
            dfxp = differentiate(differentiate(dx, 1),2)(zerosTN)

            v = dfx

            LP_Function!(F, f, dfx, xf[i], pf[i, :], v, w, t)
            LP_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w)

            j = 1

            while j <= newtonite && norm(F) > newtontol

                Δ .= - inv(J)*F

                xf[i] += Δ[1]
                pf[i, indice] += Δ[2]
                v += Δ[3]

                dx = f(xf[i] + S[1], pf[i, :] .+ Sp, t + r)

                dfx = differentiate(dx, 1)(zerosTN)
                dfp = differentiate(dx, 2)(zerosTN)
                dfxx = differentiate(differentiate(dx, 1),1)(zerosTN)
                dfxp = differentiate(differentiate(dx, 1),2)(zerosTN)

                LP_Function!(F, f, dfx, xf[i], pf[i, :], v, w, t)
                LP_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w)

                j += 1

                println("$(j) : $(norm(F))")

            end

        end

    end
    
end

#-

function Limit_Points!(f!::Function, pf::Vector{Float64}, xf::Matrix{Float64}, test::Vector{Int64}; newtonite = 8, newtontol = 1.e-16)

    m, n = size(xf)

    F = zeros(2*n + 1)
    J = zeros(2*n + 1, 2*n + 1)

    S = set_variables("s", numvars = n + 1, order = 2)

    r = zero(S[1])

    dx = [S[1] for i in 1:n]

    dfx = zeros(n,n)
    dfp = zeros(n)
    dfxx = zeros(n,n,n)
    dfxp = zeros(n,n)

    v = zeros(n)

    w = zeros(n)

    zerosTN = zeros(n + 1)

    Δ = zeros(2*n + 1)

    if m != 0

        for i in 1:m

            # println(" \n Fold $(i) : \n ")

            f!(dx, xf[i, :] .+ S[1:n], pf[i] + S[n + 1], t + r)

            dfx  .= [differentiate(dx[i], j)(zerosTN) for i in 1:n, j in 1:n]
            dfp  .= [differentiate(dx[i], n + 1)(zerosTN) for i in 1:n]
            dfxx .= [differentiate(differentiate(dx[i], j), k)(zerosTN) for i in 1:n, j in 1:n, k in 1:n]
            dfxp .= [differentiate(differentiate(dx[i], j), n + 1)(zerosTN) for i in 1:n, j in 1:n]

            v .= real.(eigvecs(dfx)[:, test[i]])

            w .= v

            LP_Function!(F, dfx, dx(zerosTN), v, w, n)
            LP_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w, n)

            j = 1

            while j < newtonite && norm(F) > newtontol

                Δ .= - inv(J)*F

                xf[i, :] .+= Δ[1:n]
                pf[i]     += Δ[n + 1]
                v        .+= Δ[n + 2:end]

                f!(dx, xf[i, :] .+ S[1:n], pf[i] + S[n + 1], t + r)

                dfx  .= [differentiate(dx[i], j)(zerosTN) for i in 1:n, j in 1:n]
                dfp  .= [differentiate(dx[i], n + 1)(zerosTN) for i in 1:n]
                dfxx .= [differentiate(differentiate(dx[i], j), k)(zerosTN) for i in 1:n, j in 1:n, k in 1:n]
                dfxp .= [differentiate(differentiate(dx[i], j), n + 1)(zerosTN) for i in 1:n, j in 1:n]

                LP_Function!(F, dfx, dx(zerosTN), v, w, n)
                LP_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w, n)

                j += 1

                # println("$(j) : $(norm(F))")

            end

        end

    end

end

#-

function Limit_Points!(f!::Function, pf::Matrix{Float64}, xf::Matrix{Float64}, test::Vector{Int64}, indice::Int64; newtonite = 8, newtontol = 1.e-16)

    m, n = size(xf)
    mp, np = size(pf)

    F = zeros(2*n + 1)
    J = zeros(2*n + 1, 2*n + 1)

    S = set_variables("s", numvars = n + 1, order = 2)
    r = zero(S[1])

    Sp = [i == indice ? S[end] : r for i in 1:np]

    dx = [S[1] for i in 1:n]

    dfx = zeros(n,n)
    dfp = zeros(n)
    dfxx = zeros(n,n,n)
    dfxp = zeros(n,n)

    zerosTN = zeros(n + 1)

    Δ = zeros(2*n + 1)

    if m != 0

        for i in 1:m

            # println(" \n Fold $(i) : \n ")

            f!(dx, xf[i, :] .+ S[1:n], pf[i, :] .+ Sp, t + r)

            dfx .= [differentiate(dx[i], j)(zerosTN) for i in 1:n, j in 1:n]
            dfp .= [differentiate(dx[i], n + 1)(zerosTN) for i in 1:n]
            dfxx .= [differentiate(differentiate(dx[i], j), k)(zerosTN) for i in 1:n, j in 1:n, k in 1:n]
            dfxp .= [differentiate(differentiate(dx[i], j), n + 1)(zerosTN) for i in 1:n, j in 1:n]

            v = real.(eigvecs(dfx)[:, test[i]])

            w = v

            LP_Function!(F, dfx, dx(zerosTN), v, w, n)
            LP_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w, n)

            j = 1

            while j < newtonite && norm(F) > newtontol

                Δ .= - inv(J)*F

                xf[i, :] .+= Δ[1:n]
                pf[i, indice] += Δ[n + 1]
                v .+= Δ[n + 2:end]

                f!(dx, xf[i, :] .+ S[1:n], pf[i, :] .+ Sp, t + r)

                dfx  .= [differentiate(dx[i], j)(zerosTN) for i in 1:n, j in 1:n]
                dfp  .= [differentiate(dx[i], n + 1)(zerosTN) for i in 1:n]
                dfxx .= [differentiate(differentiate(dx[i], j), k)(zerosTN) for i in 1:n, j in 1:n, k in 1:n]
                dfxp .= [differentiate(differentiate(dx[i], j), n + 1)(zerosTN) for i in 1:n, j in 1:n]

                LP_Function!(F, dfx, dx(zerosTN), v, w, n)
                LP_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w, n)

                j += 1

                # println("$(j) : $(norm(F))")

            end

        end

    end

end

#-

function Hopf_Points!(f!::Function, pb::Vector{Float64}, xb::Matrix{Float64}, test::Vector{Int64}; newtonite = 8, newtontol = 1.e-16)

    m, n = size(xb)

    F = zeros(2*n + 2)
    J = zeros(2*n + 2, 2*n + 2)

    S = set_variables("s", numvars = n + 1, order = 2)

    r = zero(S[1])

    dx = [S[1] for i in 1:n]

    dfx = zeros(n,n)
    dfp = zeros(n)
    dfxx = zeros(n,n,n)
    dfxp = zeros(n,n)

    zerosTN = zeros(n + 1)

    v = ones(n)

    w = zeros(n)

    Δ = zeros(2*n + 2)

    if m != 0

        for i in 1:m

            # println(" \n Hopf $(i) : \n ")

            f!(dx, xb[i, :] .+ S[1:n], pb[i] + S[n + 1], t + r)

            dfx .= [differentiate(dx[i], j)(zerosTN) for i in 1:n, j in 1:n]
            dfp .= [differentiate(dx[i], n + 1)(zerosTN) for i in 1:n]
            dfxx .= [differentiate(differentiate(dx[i], j), k)(zerosTN) for i in 1:n, j in 1:n, k in 1:n]
            dfxp .= [differentiate(differentiate(dx[i], j), n + 1)(zerosTN) for i in 1:n, j in 1:n]

            k = imag(eigvals(dfx)[test[m]])^2

            for i in 1:10
                v .= (dfx^2 + k*I(n))\v
                v .= v/norm(v)
            end

            w .= nullspace(v')[:,1]

            Hopf_Function!(F, dx(zerosTN), dfx, v, w, k, n)
            Hopf_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w, k, n)

            j = 1

            while j <= newtonite && norm(F) > newtontol

                Δ .= - inv(J)*F

                xb[i, :] .+= Δ[1:n]
                pb[i] += Δ[n + 1]
                v .+= Δ[n + 2:2*n + 1]
                k += Δ[end]

                f!(dx, xb[i, :] .+ S[1:n], pb[i] + S[n + 1], t + r)

                dfx .= [differentiate(dx[i], j)(zerosTN) for i in 1:n, j in 1:n]
                dfp .= [differentiate(dx[i], n + 1)(zerosTN) for i in 1:n]
                dfxx .= [differentiate(differentiate(dx[i], j), k)(zerosTN) for i in 1:n, j in 1:n, k in 1:n]
                dfxp .= [differentiate(differentiate(dx[i], j), n + 1)(zerosTN) for i in 1:n, j in 1:n]

                Hopf_Function!(F, dx(zerosTN), dfx, v, w, k, n)
                Hopf_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w, k, n)

                j += 1

                # println("$(j) : $(norm(F))")

            end

        end

    end

end

#-


function Hopf_Points!(f!::Function, pb::Matrix{Float64}, xb::Matrix{Float64}, test::Vector{Int64}, indice::Int64; newtonite = 8, newtontol = 1.e-16)

    m, n = size(xb)
    mp, np = size(pb)

    F = zeros(2*n + 2)
    J = zeros(2*n + 2, 2*n + 2)

    S = set_variables("s", numvars = n + 1, order = 2)
    r = zero(S[1])

    Sp = [i == indice ? S[end] : r for i in 1:np]

    dx = [S[1] for i in 1:n]

    dfx = zeros(n,n)
    dfp = zeros(n)
    dfxx = zeros(n,n,n)
    dfxp = zeros(n,n)

    zerosTN = zeros(n + 1)

    Δ = zeros(2*n + 2)

    if m != 0

        for i in 1:m

            # println(" \n Hopf $(i) : \n ")

            f!(dx, xb[i, :] .+ S[1:n], pb[i, :] .+ Sp, t + r)

            dfx .= [differentiate(dx[i], j)(zerosTN) for i in 1:n, j in 1:n]
            dfp .= [differentiate(dx[i], n + 1)(zerosTN) for i in 1:n]
            dfxx .= [differentiate(differentiate(dx[i], j), k)(zerosTN) for i in 1:n, j in 1:n, k in 1:n]
            dfxp .= [differentiate(differentiate(dx[i], j), n + 1)(zerosTN) for i in 1:n, j in 1:n]

            v = one.(xb[1, :])

            k = imag(eigvals(dfx)[test[i]])^2

            for i in 1:10
                v .= (dfx^2 + k*I(n))\v
                v .= v/norm(v)
            end

            w = nullspace(v')[:,1]

            Hopf_Function!(F, dx(zerosTN), dfx, v, w, k, n)
            Hopf_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w, k, n)

            j = 1

            while j <= newtonite && norm(F) > newtontol

                Δ .= - inv(J)*F

                xb[i, :] .+= Δ[1:n]
                pb[i, indice] += Δ[n + 1]
                v .+= Δ[n + 2:2*n + 1]
                k += Δ[end]

                f!(dx, xb[i, :] .+ S[1:n], pb[i, :] .+ Sp, t + r)

                dfx .= [differentiate(dx[i], j)(zerosTN) for i in 1:n, j in 1:n]
                dfp .= [differentiate(dx[i], n + 1)(zerosTN) for i in 1:n]
                dfxx .= [differentiate(differentiate(dx[i], j), k)(zerosTN) for i in 1:n, j in 1:n, k in 1:n]
                dfxp .= [differentiate(differentiate(dx[i], j), n + 1)(zerosTN) for i in 1:n, j in 1:n]

                Hopf_Function!(F, dx(zerosTN), dfx, v, w, k, n)
                Hopf_Jacobian!(J, dfx, dfp, dfxx, dfxp, v, w, k, n)

                j += 1

                # println("$(j) : $(norm(F))")

            end

        end

    end

end

#-


# ### Referencias

#- 

# - http://www.maths.liv.ac.uk/~bnvasiev/Past%20students/Caitlin_399.pdf
# 
# - Doedel, E.J. (2007). Lecture Notes on Numerical Analysis of Nonlinear Equations.
#   In: Krauskopf, B., Osinga, H.M., Galán-Vioque, J. (eds) Numerical Continuation 
#   Methods for Dynamical Systems. Understanding Complex Systems. Springer, Dordrecht. 
#   https://doi.org/10.1007/978-1-4020-6356-5_1 