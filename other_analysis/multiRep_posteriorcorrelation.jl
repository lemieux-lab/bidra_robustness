using DataFrames, HDF5, JLD2
include("../utils.jl")

dt = ARGS[1]
overwrite = parse(Bool, ARGS[2])

data_prefix = "public_datasets/"
bidra_params = ["LDR", "HDR", "ic50", "slope", "aac"]

function doCorrelation(df::DataFrame, rep::Int)
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
        medianCorr_df[:, :rep] = [rep]
        
        ### QQ posterior correlation
        prQQ = combine(groupby(df, :exp_id_rep1), pr_rep1 => sort => :sorted_1, pr_rep2 => sort => :sorted_2)
        
        qqCorr_df = correlationAnalysis(prQQ.sorted_1, prQQ.sorted_2)
        qqCorr_df[:, :N] = [length(unique(df.exp_id_rep1))]
        qqCorr_df[:, :dataset] = [dt]
        qqCorr_df[:, :param] = [pr]
        qqCorr_df[:, :rep] = [rep]

        if rep == 1 && overwrite
            CSV.write(results_prefix*"multiRep_medianCorrelations.csv", medianCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","rep"])
            CSV.write(results_prefix*"multiRep_qqCorrelations.csv", qqCorr_df, delim=",", appen=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","rep"])
        else
            CSV.write(results_prefix*"multiRep_medianCorrelations.csv", medianCorr_df, delim=",", append=true)
            CSV.write(results_prefix*"multiRep_qqCorrelations.csv", qqCorr_df, delim=",", append=true)
        end
    end
end

grouped_expID = load(joinpath(data_prefix, "rep_more2_pairing.jld2"))[dt]
R=10000

for r in 1:R
    println("Resampling $r")
    tmp = mapreduce(g -> DataFrame(g[rand(1:nrow(g), 2), :]), vcat, groupby(grouped_expID, :pairs))
    pairing_df = DataFrame(rep_1 = tmp.exp_id_unique[1:2:end], rep_2 = tmp.exp_id_unique[2:2:end])

    posterior_rep1 = getPosterior_h5(dt, false, Array(pairing_df.rep_1))
    posterior_rep2 = getPosterior_h5(dt, false, Array(pairing_df.rep_2))

    rename!(posterior_rep1, map(x -> "$x"*"_rep1", names(posterior_rep1)))
    rename!(posterior_rep2, map(x -> "$x"*"_rep2", names(posterior_rep2)))

    posteriorPaired_df = hcat(posterior_rep1, posterior_rep2)
    @time doCorrelation(posteriorPaired_df, r)
end