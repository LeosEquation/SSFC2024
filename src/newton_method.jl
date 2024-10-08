

function NewtonPeriodic!(i::Int64, f!::Function, F::Vector{Float64}, J::Matrix{Float64}, Δ::Vector{Float64}, x::Matrix{Float64}, p::Vector{Float64}, T::Vector{Float64}, xs::Vector{Float64}, ps::Float64, Ts::Float64, Δs::Float64, dx0::Vector{Float64}, n::Int64, dx::Vector{Float64}, time::Vector{Float64}, sol::Matrix{TaylorN{Float64}}, Sx::Vector{TaylorN{Float64}}, Sp::TaylorN{Float64}, zerosTN::Vector{Float64}, newtonite::Int64, newtontol::Float64, integorder::Int64, integtol::Float64, integmaxsteps::Int64)

    for j in 1:n 
        x[i,j] = x[i-1,j] + Δs*xs[j]
    end

    p[i] = p[i-1] + Δs*ps
    T[i] = T[i-1] + Δs*Ts

    time[2] = T[i]

    sol .= taylorinteg(f!, x[i,:] .+ Sx, time, integorder, integtol, p[i] + Sp; maxsteps = integmaxsteps, parse_eqs = false)

    Function_Periodic_Orbits!(F, x[i, :], p[i], T[i], x[i-1, :], p[i-1], T[i-1], xs, ps, Ts, Δs, sol[2, :](zerosTN), dx0, n)
    Jacobian_Periodic_Orbits!(J, f!, p[i], T[i], xs, ps, Ts, sol[2,:], dx, dx0, zerosTN, n)

    j = 1

    while j < newtonite && norm(F) > newtontol
            
        Δ .= J \ F

        for k in 1:n 
            x[i,  k] -= Δ[k]
        end

        p[i] -= Δ[n+1]
        T[i] -= Δ[n+2]

        time[2] = T[i]

        sol .= taylorinteg(f!, x[i,:] .+ Sx, time, integorder, integtol, p[i] .+ Sp; maxsteps = integmaxsteps, parse_eqs = false)

        Function_Periodic_Orbits!(F, x[i, :], p[i], T[i], x[i-1, :], p[i-1], T[i-1], xs, ps, Ts, Δs, sol[2, :](zerosTN), dx0, n)
        Jacobian_Periodic_Orbits!(J, f!, p[i], T[i], xs, ps, Ts, sol[2, :], dx, dx0, zerosTN, n)
        
        j += 1

    end

end


function NewtonPeriodic!(i::Int64, f!::Function, F::Vector{Float64}, J::Matrix{Float64}, Δ::Vector{Float64}, x::Matrix{Float64}, p::Matrix{Float64}, T::Vector{Float64}, xs::Vector{Float64}, ps::Float64, Ts::Float64, Δs::Float64, dx0::Vector{Float64}, indice::Int64, n::Int64, dx::Vector{Float64}, time::Vector{Float64}, sol::Matrix{TaylorN{Float64}}, Sx::Vector{TaylorN{Float64}}, Sp::Vector{TaylorN{Float64}}, zerosTN::Vector{Float64}, newtonite::Int64, newtontol::Float64, integorder::Int64, integtol::Float64, integmaxsteps::Int64)

    for j in 1:n 
        x[i,j] = x[i-1,j] + Δs*xs[j]
    end

    p[i, indice] = p[i-1, indice] + Δs*ps
    T[i] = T[i-1] + Δs*Ts

    time[2] = T[i]

    sol .= taylorinteg(f!, x[i,:] .+ Sx, time, integorder, integtol, p[i, :] .+ Sp; maxsteps = integmaxsteps, parse_eqs = false)

    Function_Periodic_Orbits!(F, x[i, :], p[i, indice], T[i], x[i-1, :], p[i-1,indice], T[i-1], xs, ps, Ts, Δs, sol[2, :](zerosTN), dx0, n)
    Jacobian_Periodic_Orbits!(J, f!, p[i, :], T[i], xs, ps, Ts, sol[2,:], dx, dx0, zerosTN, n)

    j = 1

    while j < newtonite && norm(F) > newtontol
            
        Δ .= J \ F

        for k in 1:n 
            x[i,  k] -= Δ[k]
        end

        p[i, indice] -= Δ[n+1]
        T[i]         -= Δ[n+2]

        time[2] = T[i]

        sol .= taylorinteg(f!, x[i,:] .+ Sx, time, integorder, integtol, p[i, :] .+ Sp; maxsteps = integmaxsteps, parse_eqs = false)

        Function_Periodic_Orbits!(F, x[i, :], p[i, indice], T[i], x[i-1, :], p[i-1,indice], T[i-1], xs, ps, Ts, Δs, sol[2, :](zerosTN), dx0, n)
        Jacobian_Periodic_Orbits!(J, f!, p[i, :], T[i], xs, ps, Ts, sol[2, :], dx, dx0, zerosTN, n)
        
        j += 1

    end

end