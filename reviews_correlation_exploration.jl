using DataFrames
# using Distributions, Statistics, StatsBase
# using Random
# using Gadfly, StatsPlots
# using Cairo, Fontconfig

include("utils.jl")

dt = "ctrpv2"#"gCSI"#"gray"
global overwrite_lm = false
global overwrite_qq = false
results_prefix = "_generated_data_reviews/correlation_reorder_pairings"
bidra_params = ["LDR", "HDR", "ic50", "slope"];

println("Dataset: ", dt)
expId_list = getExpId_h5(dt);
print("Number of experiments: ", length(expId_list), "\n")

pairings_df = getPairings_h5(dt)
println("Number of pairs: ", nrow(pairings_df), "\n")

### Make new pairings df based 
function random_pairing_order(pair)
    return sample(Array(pair), 2, replace=false)
end

### Calculate correlations
function do_correlation_lsqfit(df::DataFrame, repetition::Int, description::String)
    for pr in bidra_params
        rep1 = Symbol(pr,"_rep1")
        rep2 = Symbol(pr,"_rep2")

        tmp = filter(rep1 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)
        tmp = filter(rep2 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), tmp)

        mlCorr_df = correlationAnalysis(tmp[:,rep1], tmp[:,rep2])
        mlCorr_df[:, :N] = [nrow(tmp)]
        mlCorr_df[:, :dataset] = [dt]
        mlCorr_df[:, :param] = [pr]
        mlCorr_df[:, :repetition] = [repetition]
        mlCorr_df[:, :description] = [description]

        if overwrite_lm
            CSV.write("$(results_prefix)/mlCorrelations.csv", mlCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","reorder_rep","description"])
        else
            CSV.write("$(results_prefix)/mlCorrelations.csv", mlCorr_df, delim=",", append=true)
        end

        global overwrite_lm = false
    end
end

function do_correlation_bidra(df::DataFrame, repetition::Int, description::String)
    for pr in bidra_params
        pr_rep1 = Symbol("$pr"*"_rep1")
        pr_rep2 = Symbol("$pr"*"_rep2")

        ### QQ posterior correlation
        prQQ = combine(groupby(df, :exp_id_rep1), pr_rep1 => sort => :sorted_1, pr_rep2 => sort => :sorted_2)
        
        qqCorr_df = correlationAnalysis(prQQ.sorted_1, prQQ.sorted_2)
        qqCorr_df[:, :N] = [length(unique(df.exp_id_rep1))]
        qqCorr_df[:, :dataset] = [dt]
        qqCorr_df[:, :param] = [pr]
        qqCorr_df[:, :repetition] = [repetition]
        qqCorr_df[:, :description] = [description]

        if overwrite_qq
            CSV.write("$(results_prefix)/qqCorrelations.csv", qqCorr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","dataset","param","reorder_rep","description"])
        else 
            CSV.write("$(results_prefix)/qqCorrelations.csv", qqCorr_df, delim=",", append=true)
        end

        global overwrite_qq = false
    end
end

println("Get data, SD, and group")
data_df = getRawData_h5(dt, false)
sd_df = combine(groupby(data_df, :exp_id), :Viability => std => :std_viability)
expId_complete = sd_df[sd_df.std_viability .>= 20, :exp_id]
expId_incomplete = sd_df[sd_df.std_viability .< 20, :exp_id]

column_names = [:rep_1, :rep_2];
R = 1;

# rep = 1;
for rep in 1:R
    reorder_pairs = map(random_pairing_order, eachrow(pairings_df));
    reorder_pairs_df = DataFrame(collect.(eachrow(stack(reorder_pairs))), column_names)

    pairingComplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_complete && y ∈ expId_complete, reorder_pairs_df)
    pairingIncomplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete && y ∈ expId_incomplete, reorder_pairs_df);
    pairingMixte_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete || y ∈ expId_incomplete, reorder_pairs_df);

    mlPaired_df = getMLestimates(dt, reorder_pairs_df)

    ### ML estimates correlation 
    mlComplete_df = getMLestimates(dt, pairingComplete_df)
    mlIncomplete_df = getMLestimates(dt, pairingIncomplete_df)
    mlMixte_df = getMLestimates(dt, pairingMixte_df)

    println("Correlation for all pairs")
    @time do_correlation_lsqfit(mlPaired_df, rep, "all pairs")

    println("Correlation for complete pairs")
    @time do_correlation_lsqfit(mlComplete_df, rep, "complete pairs")

    println("Correlation for incomplete pairs")
    @time do_correlation_lsqfit(mlIncomplete_df, rep, "incomplete pairs")

    println("Correlation for mixte pairs")
    @time do_correlation_lsqfit(mlMixte_df, rep, "mixte pairs")

    ### Posterior correlation
    function make_posterior_df(dt::String, pairings::DataFrame)
        posterior_rep1 = getPosterior_h5(dt, false, String.(pairings.rep_1))
        posterior_rep2 = getPosterior_h5(dt, false, String.(pairings.rep_2))

        rename!(posterior_rep1, map(x -> String(x)*"_rep1", names(posterior_rep1)))
        rename!(posterior_rep2, map(x -> String(x)*"_rep2", names(posterior_rep2)))

        paired_df = hcat(posterior_rep1, posterior_rep2)

        return paired_df
    end

    posteriorPaired_df = make_posterior_df(dt, reorder_pairs_df);
    posteriorComplete_df = make_posterior_df(dt, pairingComplete_df);
    posteriorIncomplete_df = make_posterior_df(dt, pairingIncomplete_df);
    posteriorMixte_df = make_posterior_df(dt, pairingMixte_df);

    println("Correlation for all pairs")
    @time do_correlation_bidra(posteriorPaired_df, rep, "all pairs")

    println("Correlation for complete pairs")
    @time do_correlation_bidra(posteriorComplete_df, rep, "complete pairs")

    println("Correlation for incomplete pairs")
    @time do_correlation_bidra(posteriorIncomplete_df, rep, "incomplete pairs")

    println("Correlation for mixte pairs")
    @time do_correlation_bidra(posteriorMixte_df, rep, "mixte pairs")
end

