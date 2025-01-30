using DataFrames, HDF5, JLD2
include("../utils.jl")

dt = ARGS[1]
overwrite = parse(Bool, ARGS[2])

data_prefix = "public_datasets/"
bidra_params = ["LDR", "HDR", "ic50", "slope", "aac"]

function doCorrelation(df::DataFrame, rep::Int)
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
        mlCorr_df[:, :rep] = [rep]

        if rep == 1 && overwrite
            CSV.write(results_prefix*"multiRep_mlCorrelations.csv", mlCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","rep"])
        else
            CSV.write(results_prefix*"multiRep_mlCorrelations.csv", mlCorr_df, delim=",", append=true)
        end
    end
end

grouped_expID = load(joinpath(data_prefix, "rep_more2_pairing.jld2"))[dt]
R = 10000




for r in 1:R
    println("Resampling $r")

    tmp = mapreduce(g -> DataFrame(g[rand(1:nrow(g), 2), :]), vcat, groupby(grouped_expID, :pairs))
    pairing_df = DataFrame(rep_1 = tmp.exp_id_unique[1:2:end], rep_2 = tmp.exp_id_unique[2:2:end])
    mlPaired_df = getMLestimates(dt, pairing_df)

    ### replace Inf aac by Nan
    mlPaired_df[!, :aac_rep_1] = replace(mlPaired_df.aac_rep1, Inf => NaN, -Inf => NaN)
    mlPaired_df[!, :aac_rep_2] = replace(mlPaired_df.aac_rep2, Inf => NaN, -Inf => NaN)

    ### Correlation on all pairs
    @time doCorrelation(mlPaired_df, r)
end

