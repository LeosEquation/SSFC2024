
module NewtonMethod

export Newton

using TaylorSeries

function Newton(f::Function,p::Union{Float64,Vector{Float64}},x0::Float64,orden::Int64; ite = 30, tol = 1e-16)
    x_new = x0
    i = 1
    while i <= ite && abs(f(x_new,p)) > tol
        x_old = x_new
        x = Taylor1([x_old,1],orden)
        x_new = x_old - (f(x,p)/differentiate(f(x,p)))(0.0)
        i += 1
    end
    if abs(f(x_new,p)) > tol
        throw(ArgumentError("El método no convergió en la tolerancia deseada"))
    end
    """if abs(f(x_new,p)) <= tol
        return x_new
    else
        #println("El método no convergió en la tolerancia deseada\t |f($(x_new),$(p))| = $(abs(f(x_new,p))) > $(tol)")
        #return NaN
    end"""
end

end