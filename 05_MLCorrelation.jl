using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")

results_prefix = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/"

### Define dataset to analyze
dt = ARGS[1]
#dt = "gCSI"
dt = "gray"
#dt = "ctrpv2"

pairings_df = getPairings(dt)
mlPaired_df = getMLestimates([dt], true, pairings_df)
bidra_params = ["LDR", "HDR", "ic50", "slope", "aac"]

### replace Inf aac by Nan
mlPaired_df[!, :aac_rep_1] = replace(mlPaired_df.aac_rep_1, Inf => NaN, -Inf => NaN)
mlPaired_df[!, :aac_rep_2] = replace(mlPaired_df.aac_rep_2, Inf => NaN, -Inf => NaN)

for pr in bidra_params
    rep1 = Symbol(pr,"_rep_1")
    rep2 = Symbol(pr,"_rep_2")

    tmp = filter(rep1 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), mlPaired_df)
    tmp = filter(rep2 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), tmp)

    mlCorr_df = correlationAnalysis(tmp[:,rep1], tmp[:,rep2])
    mlCorr_df[:, :N] = [nrow(tmp)]
    mlCorr_df[:, :dataset] = [dt*"_all"]
    mlCorr_df[:, :param] = [pr]

    CSV.write(results_prefix*"mlCorrelations.csv", mlCorr_df, delim=",", append=true)
end

subset_converged = filter([:convergence_rep_1, :convergence_rep_2] => (a, b) -> a + b == 2, mlPaired_df)
for pr in bidra_params
    rep1 = Symbol(pr,"_rep_1")
    rep2 = Symbol(pr,"_rep_2")

    tmp = filter(rep1 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), subset_converged)
    tmp = filter(rep2 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), tmp)

    mlCorr_df = correlationAnalysis(tmp[:,rep1], tmp[:,rep2])
    mlCorr_df[:, :N] = [nrow(tmp)]
    mlCorr_df[:, :dataset] = [dt*"_converge"]
    mlCorr_df[:, :param] = [pr]

    CSV.write(results_prefix*"mlCorrelations.csv", mlCorr_df, delim=",", append=true)
end
