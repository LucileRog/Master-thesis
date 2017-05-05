pwd()

include("src/autarkyeq.jl")
include("src/tradeeq.jl")
include("src/subsistconst.jl")

autarkyeq.runall(Lmax::Flaot64, t::Float64, sigma::Flaot64)
autarkyplots.plotall()

tradeeqdf.runall(Lmax::Vector, T::Vector, sigma::Vector)

subsistconst.subsist(Lmax, T, sigma)
