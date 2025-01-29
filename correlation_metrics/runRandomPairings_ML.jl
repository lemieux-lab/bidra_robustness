using DataFrames
using Statistics, StatsBase

include("../utils.jl")

### Gloab variables
results_prefix = "_generated_data/"
datasets = ["gray", "gCSI", "ctrpv2"]
bidra_params = ["LDR", "HDR", "ic50", "slope"]

function randomPairings(dt::String, n::Int64, expID_list::Array, rep::Int64, description::String, addHeader::Bool)
    println("------  Generate random exp_id list")
    rdn_sort = sample(1:n, n, replace = false)
    exp_id_rdm = expID_list[rdn_sort]
    rdm_pairing = DataFrame(rep_1=exp_id_rdm[1:2:end], rep_2=exp_id_rdm[2:2:end])

    df_random_paired = getMLestimates([dt], rdm_pairing)
    for pr in bidra_params
        mlCorr_df = correlationAnalysis(df_random_paired[:, "$pr"*"_rep1"], df_random_paired[:, "$pr"*"_rep2"])
        mlCorr_df[:, :N] = [n]
        mlCorr_df[:, :dataset] = [dt]
        mlCorr_df[:, :param] = [pr]
        mlCorr_df[:, :rep] = [rep]
        mlCorr_df[:, :description] = [description]
        
        if addHeader
            CSV.write(results_prefix*"mlRandomCorrelations.csv", mlCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","rep","method","description"])
        else
            CSV.write(results_prefix*"mlRandomCorrelations.csv", mlCorr_df, delim=",", append=true)
        end
    end
end

function getExpIdList(pairings::DataFrame)
    return vcat(pairings.rep_1, pairings.rep_2)
end


dt_count = 1
for dt in datasets
    println("### $dt")

    println("Get data")
    data_df = getRawData_h5(dt, false)
    sd_df = combine(groupby(data_df, :exp_id), :Viability => std => :std_viability)
    expId_complete = sd_df[sd_df.std_viability .>= 20, :exp_id]
    expId_incomplete = sd_df[sd_df.std_viability .< 20, :exp_id]

    println("Get pairings")
    pairings_df = getPairings_h5(dt, si)
    pairingComplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_complete && y ∈ expId_complete, pairings_df)
    pairingIncomplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete && y ∈ expId_incomplete, pairings_df)
    pairingMixte_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete || y ∈ expId_incomplete, pairings_df)

    println("Get ML estimates")
    mlPaired_df = getMLestimates([dt], pairings_df)
    pairingConverge_df = filter([:convergence_rep1, :convergence_rep2] => (a, b) -> a + b == 2, mlPaired_df)

    N = [nrow(df)*2 for df in [pairings_df, pairingComplete_df, pairingIncomplete_df, pairingMixte_df, pairingConverge_df]]

    expId_list_all = getExpIdList(pairings_df)
    expId_list_complete = getExpIdList(pairingComplete_df)
    expId_list_incomplete = getExpIdList(pairingIncomplete_df)
    expId_list_mixte = getExpIdList(pairingMixte_df)
    expId_list_converge = getExpIdList(pairingConverge_df)

    ### Radomly paired correlation analysis
    for i in 1:100
        println("Replicated sampling #", i)
        println("--- All pairs")  

        if dt_count + i == 1
            randomPairings(dt, N[1], expId_list_all, i, "all pairs", true)
        else
            randomPairings(dt, N[1], expId_list_all, i, "all pairs", false)
        end

        println("--- Complete pairs")
        randomPairings(dt, N[2], expId_list_complete, i, "complete pairs", false)
        println("--- Incomplete pairs")
        randomPairings(dt, N[3], expId_list_incomplete, i, "incomplete pairs", false)
        println("--- Mixte pairs")
        randomPairings(dt, N[4], expId_list_mixte, i, "mixte pairs", false)
        println("--- Converge pairs")
        randomPairings(dt, N[5], expId_list_converge, i, "converge pairs", false)
    end

    dt_count += 1
end
