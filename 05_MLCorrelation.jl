using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")

results_prefix = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/"

gCSIpairings_df = getPairings("gCSI")
gCSImlPaired_df = getMLestimates(["gCSI"], true, gCSIpairings_df)
bidra_params = ["LDR", "HDR", "ic50", "slope", "aac"]

### replace Inf aac by Nan
gCSImlPaired_df[!, :aac_rep_1] = replace(gCSImlPaired_df.aac_rep_1, Inf => NaN)
gCSImlPaired_df[!, :aac_rep_2] = replace(gCSImlPaired_df.aac_rep_2, Inf => NaN)

for pr in bidra_params
    rep1 = Symbol(pr,"_rep_1")
    rep2 = Symbol(pr,"_rep_2")

    tmp = filter(rep1 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), gCSImlPaired_df)
    tmp = filter(rep2 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), tmp)

    mlCorr_df = correlationAnalysis(tmp[:,rep1], tmp[:,rep2])
    mlCorr_df[:, :N] = [nrow(tmp)]
    mlCorr_df[:, :dataset] = ["gCSI_all"]
    mlCorr_df[:, :param] = [pr]

    CSV.write(results_prefix*"mlCorrelations.csv", mlCorr_df, delim=",", append=true)
end

gCSIsubset_converged = filter([:convergence_rep_1, :convergence_rep_2] => (a, b) -> a + b == 2, gCSImlPaired_df)
for pr in bidra_params
    rep1 = Symbol(pr,"_rep_1")
    rep2 = Symbol(pr,"_rep_2")

    tmp = filter(rep1 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), gCSIsubset_converged)
    tmp = filter(rep2 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), tmp)

    mlCorr_df = correlationAnalysis(tmp[:,rep1], tmp[:,rep2])
    mlCorr_df[:, :N] = [nrow(tmp)]
    mlCorr_df[:, :dataset] = ["gCSI_converge"]
    mlCorr_df[:, :param] = [pr]

    CSV.write(results_prefix*"mlCorrelations.csv", mlCorr_df, delim=",", append=true)
end
