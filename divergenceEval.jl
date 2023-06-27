using Turing, MCMCChains
using JLD2, FileIO
using Serialization

include("utils.jl")

for dt in ["gray", "gCSI", "ctrpv2"]
    dt_goodness_fit = get_divergenceRate(dt)

    println("$(dt) exp. with divergence greater than 0: $(sum(dt_goodness_fit.divergence .!= 0.0))")
end