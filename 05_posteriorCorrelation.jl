using DataFrames
using Statistics, StatsBase

include("utils.jl")

results_prefix = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/"

gCSIpairings_df = getPairings("gCSI")
gCSIposteriorPaired_df = getPairedPosterior("gCSI", gCSIpairings_df)
bidra_params = ["LDR", "HDR", "ic50", "slope"]

for pr in bidra_params
    prSubset = filter(:param => x -> x == pr, gCSIposteriorPaired_df)

    ### Median correlation
    prMedian = combine(groupby(prSubset, :exp_id1), :exp1_val => median => :median_1, :exp2_val => median => :median_2)
    
    medianCorr_df = correlationAnalysis(prMedian.median_1, prMedian.median_2)
    medianCorr_df[:, :N] = [nrow(prMedian)]
    medianCorr_df[:, :dataset] = ["gCSI"]
    medianCorr_df[:, :param] = [pr]
    
    CSV.write(results_prefix*"medianCorrelations.csv", medianCorr_df, delim=",", append=true)


    ### Posterior correlation
    posteriorCorr_df = correlationAnalysis(prSubset.exp1_val, prSubset.exp2_val)
    posteriorCorr_df[:, :N] = [length(unique(prSubset.exp_id1))]
    posteriorCorr_df[:, :dataset] = ["gCSI"]
    posteriorCorr_df[:, :param] = [pr]

    CSV.write(results_prefix*"posteriorCorrelations.csv", posteriorCorr_df, delim=",", append=true)


    ### CI posterior correlation
    quantilesCI = [2.5, 50.0, 97.5]
    prCI = combine(groupby(prSubset, :exp_id1), :exp1_val => (x -> percentile(x, quantilesCI)) => :quantile_1, :exp2_val => (x -> percentile(x, quantilesCI)) => :quantile_2)

    ciCorr_df = correlationAnalysis(prCI.quantile_1, prCI.quantile_2)
    ciCorr_df[:, :N] = [length(unique(prSubset.exp_id1))]
    ciCorr_df[:, :dataset] = ["gCSI"]
    ciCorr_df[:, :param] = [pr]

    CSV.write(results_prefix*"ciCorrelations.csv", ciCorr_df, delim=",", append=true)

    ### QQ posterior correlation
    prQQ = combine(groupby(prSubset, :exp_id1), :exp1_val => sort => :sorted_1, :exp2_val => sort => :sorted_2)
    
    qqCorr_df = correlationAnalysis(prQQ.sorted_1, prQQ.sorted_2)
    qqCorr_df[:, :N] = [length(unique(prSubset.exp_id1))]
    qqCorr_df[:, :dataset] = ["gCSI"]
    qqCorr_df[:, :param] = [pr]

    CSV.write(results_prefix*"qqCorrelations.csv", qqCorr_df, delim=",", append=true)
end