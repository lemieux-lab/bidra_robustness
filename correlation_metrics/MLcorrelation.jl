using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("../utils.jl")

### Define dataset to analyze
dt = ARGS[1]
overwrite = parse(Bool, ARGS[2])

@time expId_list = getExpId_h5(dt);
@time si = StrIndex(expId_list);
pairings_df = getPairings_h5(dt, si)
mlPaired_df = getMLestimates(si, pairings_df)
bidra_params = ["LDR", "HDR", "ic50", "slope", "aac"]

### replace Inf aac by Nan
mlPaired_df[!, :aac_rep1] = replace(mlPaired_df.aac_rep1, Inf => NaN, -Inf => NaN)
mlPaired_df[!, :aac_rep2] = replace(mlPaired_df.aac_rep2, Inf => NaN, -Inf => NaN)

function doCorrelation(df::DataFrame, description::String)
    results_prefix = "_generated_data/"

    for pr in bidra_params
        rep1 = Symbol(pr,"_rep1")
        rep2 = Symbol(pr,"_rep2")

        tmp = filter(rep1 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)
        tmp = filter(rep2 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), tmp)

        mlCorr_df = correlationAnalysis(tmp[:,rep1], tmp[:,rep2])
        mlCorr_df[:, :N] = [nrow(tmp)]
        mlCorr_df[:, :dataset] = [dt]
        mlCorr_df[:, :param] = [pr]
        mlCorr_df[:, :description] = [description]

        if overwrite
            CSV.write(results_prefix*"mlCorrelations.csv", mlCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","description"])
        else
            CSV.write(results_prefix*"mlCorrelations.csv", mlCorr_df, delim=",", append=true)
        end
    end
end


println("Correlation for all pairs")
@time doCorrelation(mlPaired_df, "all pairs")

println("Correlation for converged pairs")
subset_converged = filter([:convergence_rep_1, :convergence_rep_2] => (a, b) -> a + b == 2, mlPaired_df)
println("There are ", nrow(subset_converged), " pairs for which both experiments converged")
@time doCorrelation(subset_converged, "converged pairs")

println()
println("Get data, SD, and group")
data_df = getRawData_h5(dt, false)
sd_df = combine(groupby(data_df, :exp_id), :Viability => std => :std_viability)
expId_complete = sd_df[sd_df.std_viability .>= 20, :exp_id]
expId_incomplete = sd_df[sd_df.std_viability .< 20, :exp_id]
pairingComplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_complete && y ∈ expId_complete, pairings_df)
pairingIncomplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete && y ∈ expId_incomplete, pairings_df)
pairingMixte_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete || y ∈ expId_incomplete, pairings_df)

mlComplete_df = getMLestimates([dt], pairingComplete_df)
mlIncomplete_df = getMLestimates([dt], pairingIncomplete_df)
mlMixte_df = getMLestimates([dt], pairingMixte_df)

println("Correlation for complete pairs")
@time doCorrelation(mlComplete_df, "complete pairs")

println("Correlation for incomplete pairs")
@time doCorrelation(mlIncomplete_df, "incomplete pairs")

println("Correlation for mixte pairs")
@time doCorrelation(mlMixte_df, "mixte pairs")