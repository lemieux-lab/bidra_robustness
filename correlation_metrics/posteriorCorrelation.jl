using DataFrames
using Statistics, StatsBase

include("../utils.jl")

### Define dataset to analyze
dt = ARGS[1]
overwrite = parse(Bool, ARGS[2])
println(overwrite)

@time expId_list = getExpId_h5(dt);
@time si = StrIndex(expId_list);
pairings_df = getPairings_h5(dt, si)
posteriorPaired_df = getPairedPosterior_h5(dt, si)
bidra_params = ["LDR", "HDR", "ic50", "slope", "aac"]

function doCorrelation(df::DataFrame, description::String)
    results_prefix = "_generated_data/"

    for pr in bidra_params
        ### Median correlation
        pr_rep1 = Symbol("$pr"*"_rep1")
        pr_rep2 = Symbol("$pr"*"_rep2")
        prMedian = combine(groupby(df, :exp_id_rep1), pr_rep1 => median => :median_1, pr_rep2 => median => :median_2)
        
        medianCorr_df = correlationAnalysis(prMedian.median_1, prMedian.median_2)
        medianCorr_df[:, :N] = [nrow(prMedian)]
        medianCorr_df[:, :dataset] = [dt]
        medianCorr_df[:, :param] = [pr]
        medianCorr_df[:, :description] = [description]
        
        if overwrite 
            CSV.write(results_prefix*"medianCorrelations.csv", medianCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","description"])
        else 
            CSV.write(results_prefix*"medianCorrelations.csv", medianCorr_df, delim=",", append=true)
        end
    
        ### Posterior correlation
        posteriorCorr_df = correlationAnalysis(df[:, pr_rep1], df[:, pr_rep2])
        posteriorCorr_df[:, :N] = [length(unique(df.exp_id_rep1))]
        posteriorCorr_df[:, :dataset] = [dt]
        posteriorCorr_df[:, :param] = [pr]
        posteriorCorr_df[:, :description] = [description]
        
        if overwrite 
            CSV.write(results_prefix*"posteriorCorrelations.csv", posteriorCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","description"])
        else 
            CSV.write(results_prefix*"posteriorCorrelations.csv", posteriorCorr_df, delim=",", append=true)
        end
    
        ### CI posterior correlation
        quantilesCI = [2.5, 50.0, 97.5]
        prCI = combine(groupby(df, :exp_id_rep1), pr_rep1 => (x -> percentile(x, quantilesCI)) => :quantile_1, pr_rep2 => (x -> percentile(x, quantilesCI)) => :quantile_2)
    
        ciCorr_df = correlationAnalysis(prCI.quantile_1, prCI.quantile_2)
        ciCorr_df[:, :N] = [length(unique(df.exp_id_rep1))]
        ciCorr_df[:, :dataset] = [dt]
        ciCorr_df[:, :param] = [pr]
        ciCorr_df[:, :description] = [description]
    
        CSV.write(results_prefix*"ciCorrelations.csv", ciCorr_df, delim=",", append=true)
    
        ### QQ posterior correlation
        prQQ = combine(groupby(df, :exp_id_rep1), pr_rep1 => sort => :sorted_1, pr_rep2 => sort => :sorted_2)
        
        qqCorr_df = correlationAnalysis(prQQ.sorted_1, prQQ.sorted_2)
        qqCorr_df[:, :N] = [length(unique(df.exp_id_rep1))]
        qqCorr_df[:, :dataset] = [dt]
        qqCorr_df[:, :param] = [pr]
        qqCorr_df[:, :description] = [description]
    
        
        if overwrite 
            CSV.write(results_prefix*"qqCorrelations.csv", qqCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","description"])
        else 
            CSV.write(results_prefix*"qqCorrelations.csv", qqCorr_df, delim=",", append=true)
        end
        
    end
end

function createPosteriorDf(dt::String, pairings::DataFrame)
    posterior_rep1 = getPosterior_h5(dt, false, si, Array(pairings.rep_1))
    posterior_rep2 = getPosterior_h5(dt, false, si, Array(pairings.rep_2))

    rename!(posterior_rep1, map(x -> "$x"*"_rep1", names(posterior_rep1)))
    rename!(posterior_rep2, map(x -> "$x"*"_rep2", names(posterior_rep2)))

    paired_df = hcat(posterior_rep1, posterior_rep2)

    return paired_df
end

println("Correlation for all pairs")
@time doCorrelation(posteriorPaired_df, "all pairs")

println()
println("Get data, SD, and group")
data_df = getRawData_h5(dt, false, si)
sd_df = combine(groupby(data_df, :exp_id), :Viability => std => :std_viability)
expId_complete = sd_df[sd_df.std_viability .>= 20, :exp_id]
expId_incomplete = sd_df[sd_df.std_viability .< 20, :exp_id]
println("--- There are ", length(expId_complete), " exp. with SD ≥ 20")
println("--- There are ", length(expId_incomplete), " exp. with SD < 20")

pairingComplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_complete && y ∈ expId_complete, pairings_df)
println("--- There are ", nrow(pairingComplete_df), " pairs with both exp. having SD ≥ 20")
posteriorComplete_df = createPosteriorDf(dt, pairingComplete_df)

pairingIncomplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete && y ∈ expId_incomplete, pairings_df)
println("--- There are ", nrow(pairingIncomplete_df), " pairs with both exp. having SD < 20")
posteriorIncomplete_df = createPosteriorDf(dt, pairingIncomplete_df)

pairingMixte_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete || y ∈ expId_incomplete, pairings_df)
println("--- There are ", nrow(pairingMixte_df), " pairs with at least one exp. having SD < 20")
posteriorMixte_df = createPosteriorDf(dt, pairingMixte_df)

println("Correlation for complete pairs")
@time doCorrelation(posteriorComplete_df, "complete pairs")

println("Correlation for incomplete pairs")
@time doCorrelation(posteriorIncomplete_df, "incomplete pairs")

println("Correlation for mixte pairs")
@time doCorrelation(posteriorMixte_df, "mixte pairs")