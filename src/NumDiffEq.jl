module NumDiffEq

export Equilibrium, Stability_intervals, Limit_Points, Hopf_function, Hopf_initial_values, Hopf_Jacobian

using TaylorSeries, LinearAlgebra

include("diff_tools.jl")
include("newton.jl")
include("implicit_function.jl")
include("psarc.jl")
include("equilibrium.jl")
include("stability.jl")
include("bifurcation.jl")

end