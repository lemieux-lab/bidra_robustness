using DataFrames
using Statistics, StatsBase

include("utils.jl")

results_prefix = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/"

dt = "gCSI"
pairings_df = getPairings(dt)
mlPaired_df = getMLestimates([dt], true, pairings_df)
bidra_params = ["LDR", "HDR", "ic50", "slope", "aac"]

mlPaired_df[!, :aac_rep_1] = replace(mlPaired_df.aac_rep_1, Inf => NaN)
mlPaired_df[!, :aac_rep_2] = replace(mlPaired_df.aac_rep_2, Inf => NaN)
subset_converged = filter([:convergence_rep_1, :convergence_rep_2] => (a, b) -> a + b == 2, mlPaired_df)

df_list = [mlPaired_df,subset_converged]
lb_list = ["all", "converged"]

#### Radomly paired correlation analysis
for i in 1:100
    println("Replicated sampling #", i)

    for D in [1, 2]
        df = df_list[D]
        exp_id = copy(df.rep_1)
        append!(exp_id, df.rep_2)

        rdn_sort = sample(1:length(exp_id), length(exp_id), replace = false)
        exp_id_rdm = exp_id[rdn_sort]
        rdm_pairing = DataFrame(rep_1=exp_id_rdm[1:2:end], rep_2=exp_id_rdm[2:2:end])

        df_mle_paired_rdm = getMLestimates([dt], true, rdm_pairing)
        df_mle_paired_rdm[!, :aac_rep_1] = replace(df_mle_paired_rdm.aac_rep_1, Inf => NaN)
        df_mle_paired_rdm[!, :aac_rep_2] = replace(df_mle_paired_rdm.aac_rep_2, Inf => NaN)

        for pr in bidra_params
            rep1 = Symbol(pr,"_rep_1")
            rep2 = Symbol(pr,"_rep_2")
        
            tmp = filter(rep1 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df_mle_paired_rdm)
            tmp = filter(rep2 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), tmp)
        
            mlCorr_df = correlationAnalysis(tmp[:,rep1], tmp[:,rep2])
            mlCorr_df[:, :N] = [nrow(tmp)]
            mlCorr_df[:, :dataset] = ["gCSI_"*lb_list[D]]
            mlCorr_df[:, :param] = [pr]
            mlCorr_df[:, :rep] = [i]
        
            CSV.write(results_prefix*"mlRandomCorrelations.csv", mlCorr_df, delim=",", append=true)
        end
    end
end
