using Turing, MCMCChains
using JLD2, FileIO
using Serialization

include("utils.jl")

for dt in ["gray", "gCSI", "ctrpv2"]
    dt_goodness_fit = get_divergenceRate(dt)
    #dt_crfs = get_csrf(dt)

    println("$(dt) exp. with divergence greater than 0: $(sum(dt_goodness_fit.divergence .!= 0.0))")

    #good_crfs = nrow(filter(:val => x -> 0.8 ≤ x ≤ 1.2, dt_crfs))

    #println("$(dt) exp. with good CSRF: $(good_crfs) ($(good_crfs/nrow(dt_crfs)))")
end