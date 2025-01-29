using DataFrames
using Statistics, StatsBase

include("../utils.jl")

### Gloab variables
results_prefix = "_generated_data/"
datasets = ["gray", "gCSI", "ctrpv2"]
bidra_params = ["LDR", "HDR", "ic50", "slope", "aac"]

function randomPairings(dt::String, n::Int64, expId_list::Array, rep::Int64, description::String, addHeader::Bool)
    println("------  Generate random exp_id list")
    rdn_sort = sample(1:n, n, replace = false)
    exp_id_rdm = expId_list[rdn_sort]
    rdm_pairing = DataFrame(exp_id1=exp_id_rdm[1:2:end], exp_id2=exp_id_rdm[2:2:end])

    println("------ Generate random pairings")
    posterior_rep1 = getPosterior_h5(dt, false, Array(rdm_pairing.exp_id1))
    rename!(posterior_rep1, map(x -> "$x"*"_rep1", names(posterior_rep1)))
    posterior_rep2 = getPosterior_h5(dt, false, Array(rdm_pairing.exp_id2))
    rename!(posterior_rep2, map(x -> "$x"*"_rep2", names(posterior_rep2)))

    df_random_paired = hcat(posterior_rep1, posterior_rep2)
    for pr in bidra_params
        println("------ ", pr)
        posteriorCorr_df = correlationAnalysis(df_random_paired[:, "$pr"*"_rep1"], df_random_paired[:, "$pr"*"_rep2"])
        posteriorCorr_df[:, :N] = [n]
        posteriorCorr_df[:, :dataset] = [dt]
        posteriorCorr_df[:, :param] = [pr]
        posteriorCorr_df[:, :rep] = [rep]
        posteriorCorr_df[:, :method] = ["posterior"]
        posteriorCorr_df[:, :description] = [description]

        if addHeader 
            CSV.write(results_prefix*"bidraRandomCorrelation.csv", posteriorCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","rep","method","description"])
        else
            CSV.write(results_prefix*"bidraRandomCorrelation.csv", posteriorCorr_df, delim=",", append=true)
        end

        prQQ = combine(groupby(df_random_paired, :exp_id_rep1), Symbol("$pr"*"_rep1") => sort => :sorted_1, Symbol("$pr"*"_rep2") => sort => :sorted_2)
        qqCorr_df = correlationAnalysis(prQQ.sorted_1, prQQ.sorted_2)
        qqCorr_df[:, :N] = [n]
        qqCorr_df[:, :dataset] = [dt]
        qqCorr_df[:, :param] = [pr]
        qqCorr_df[:, :rep] = [rep]
        qqCorr_df[:, :method] = ["qq"]
        qqCorr_df[:, :description] = [description]
        CSV.write(results_prefix*"bidraRandomCorrelation.csv", qqCorr_df, delim=",", append=true)
    end
end

function getExpIdList(pairings::DataFrame)
    return vcat(pairings.rep_1, pairings.rep_2)
end

dt_count = 1
for dt in datasets
    println("### $dt")

    println("Get data")
    @time expId_list = getExpId_h5(dt);
    @time si = StrIndex(expId_list);
    data_df = getRawData_h5(dt, false, si)
    sd_df = combine(groupby(data_df, :exp_id), :Viability => std => :std_viability)
    expId_complete = sd_df[sd_df.std_viability .>= 20, :exp_id]
    expId_incomplete = sd_df[sd_df.std_viability .< 20, :exp_id]

    println("Get pairings")
    pairings_df = getPairings_h5(dt, si)
    pairingComplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_complete && y ∈ expId_complete, pairings_df)
    pairingIncomplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete && y ∈ expId_incomplete, pairings_df)
    pairingMixte_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete || y ∈ expId_incomplete, pairings_df)

    N = [nrow(df)*2 for df in [pairings_df, pairingComplete_df, pairingIncomplete_df, pairingMixte_df]]

    expId_list_all = getExpIdList(pairings_df)
    expId_list_complete = getExpIdList(pairingComplete_df)
    expId_list_incomplete = getExpIdList(pairingIncomplete_df)
    expId_list_mixte = getExpIdList(pairingMixte_df)
    
    
    ### Radomly paired correlation analysis for all pairs
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
        println("--- Mixtes pairs")
        randomPairings(dt, N[4], expId_list_mixte, i, "mixte pairs", false)
    end

    dt_count += 1
end
