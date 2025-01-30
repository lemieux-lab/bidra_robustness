using DataFrames, HDF5, JLD2
using Gadfly, StatsPlots, CairoMakie
#using JuBox

include("../utils.jl")

info_prefix = "public_datasets/curves_info/"
all_expId = mapreduce(dt -> getExpId_h5(dt), vcat, ["gray", "gCSI", "ctrpv2"])
si_expId = StrIndex(all_expId)
bidra_params = ["LDR", "HDR", "ic50", "slope"]
scaleColor = Scale.lab_gradient("gray95","black")

### Get info of all datasets
function getInfo(dt::String)
    info_df = readCSV(info_prefix*dt*"_info.csv", true)[:, 1:3]
    info_df = replaceChar(info_df, :exp_id)

    ### Slect experiments given our other analysis
    expId_list = getExpId_h5(dt)
    info_df = filter(:exp_id => x -> x ∈ expId_list, info_df)
    info_df[!, :exp_id] = [si_expId.str2id[v] for v in info_df.exp_id]

    info_df[!, "pairs"] = string.(info_df.drugid, ":", info_df.cellid)
    replicates_df = combine(groupby(info_df, :pairs), :exp_id=>(length∘unique)=>:n_rep)
    singleton_pairs = filter(:n_rep => x -> x == 1, replicates_df)[:, :pairs]

    return filter(:pairs => x -> x ∈ singleton_pairs, info_df)[:, [:exp_id, :pairs]]
end

function createPairedList(df::DataFrame, dt::Array{String, 1})
    rename!(df, Symbol("exp_id_"*dt[1]) => "rep_1")
    rename!(df, Symbol("exp_id_"*dt[2]) => "rep_2")
    return df[:, [:rep_1, :rep_2]]
end

function plotScatter(df::DataFrame, pair::String)
    figure_prefix = "_generated_figures/supp_fig/across_dataset/"
    for pr in bidra_params
        p = Gadfly.plot(df, x=Symbol(pr*"_rep1"), y=x=Symbol(pr*"_rep2"), Geom.point(), Geom.abline())
        draw(PDF(figure_prefix*pair*"_LM_scatter__$pr.pdf", 4inch, 4inch), p)
    end
end

function plotHexbin(df::DataFrame, pair::String)
    figure_prefix = "_generated_figures/supp_fig/across_dataset/"
    for pr in bidra_params
        p = Gadfly.plot(df, x=Symbol(pr*"_rep1"), y=x=Symbol(pr*"_rep2"), Geom.hexbin(xbincount=80, ybincount=80), Geom.abline(),
        Scale.color_continuous(colormap=scaleColor, minvalue=1))
        draw(PDF(figure_prefix*pair*"_ML_hexbin_$pr.pdf", 4inch, 4inch), p)
    end
end

function doCorrelationML(df::DataFrame, pair::String)
    results_prefix = "_generated_data/"

    for pr in bidra_params
        rep1 = Symbol(pr,"_rep1")
        rep2 = Symbol(pr,"_rep2")

        tmp = filter(rep1 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)
        tmp = filter(rep2 => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), tmp)

        mlCorr_df = correlationAnalysis(tmp[:,rep1], tmp[:,rep2])
        mlCorr_df[:, :N] = [nrow(tmp)]
        mlCorr_df[:, :param] = [pr]
        mlCorr_df[:, :method] = ["ML"]
        mlCorr_df[:, :pair] = [pair]

        CSV.write(results_prefix*"acrossDataset_correlations.csv", mlCorr_df, delim=",", append=true)
    end
end

function doCorrelationPosterior(df::DataFrame, pair::String)
    results_prefix = "_generated_data/"

    for pr in bidra_params
        ### Median correlation
        pr_rep1 = Symbol("$pr"*"_rep1")
        pr_rep2 = Symbol("$pr"*"_rep2")
        prMedian = combine(groupby(df, :exp_id_rep1), pr_rep1 => median => :median_1, pr_rep2 => median => :median_2)
        
        medianCorr_df = correlationAnalysis(prMedian.median_1, prMedian.median_2)
        medianCorr_df[:, :N] = [nrow(prMedian)]
        medianCorr_df[:, :param] = [pr]
        medianCorr_df[:, :method] = ["median"]
        medianCorr_df[:, :pair] = [pair]

        CSV.write(results_prefix*"acrossDataset_correlations.csv", medianCorr_df, delim=",", append=true)
    
        ### QQ posterior correlation
        prQQ = combine(groupby(df, :exp_id_rep1), pr_rep1 => sort => :sorted_1, pr_rep2 => sort => :sorted_2)
        
        qqCorr_df = correlationAnalysis(prQQ.sorted_1, prQQ.sorted_2)
        qqCorr_df[:, :N] = [length(unique(df.exp_id_rep1))]
        qqCorr_df[:, :param] = [pr]
        qqCorr_df[:, :method] = ["QQ"]
        qqCorr_df[:, :pair] = [pair]
    
        CSV.write(results_prefix*"acrossDataset_correlations.csv", qqCorr_df, delim=",", append=true)

        p = Gadfly.plot(prQQ, x=:sorted_1, y=:sorted_2, Geom.hexbin(xbincount=80, ybincount=80), Geom.abline(),
                        Scale.color_continuous(colormap=scaleColor, minvalue=1))
        draw(PDF(figure_prefix*pair*"_posterior_$pr.pdf", 4inch, 4inch), p)
    end
end

function singletonCorrelation(common_df::DataFrame, pair_lst::Array{String, 1})
    println("Assessing correlation of: ", pair_lst)
    println("--- ML")
    pairing_df = createPairedList(copy(common_df), pair_lst)
    mlPaired_df = getMLestimates(si_expId, pairing_df)

    pair = pair_lst[1]*"_"*pair_lst[2]
    doCorrelationML(mlPaired_df, pair)
    plotScatter(mlPaired_df, pair)
    plotHexbin(mlPaired_df, pair)

    println("--- BiDRA")
    posteriorPaired_df = getPairedPosterior_h5(pairing_df, si_expId, pair_lst)
    doCorrelationPosterior(posteriorPaired_df, pair)
end

### Individual singletons by dataset
singleton_gCSI = getInfo("gCSI")
singleton_gray = getInfo("gray")
singleton_ctrpv2 = getInfo("ctrpv2")

### Common singletons by pairs of datasets
common_gCSI_gray = innerjoin(singleton_gCSI, singleton_gray, on=:pairs, renamecols="_gCSI" => "_gray")
common_gCSI_ctrpv2 = innerjoin(singleton_gCSI, singleton_ctrpv2, on=:pairs, renamecols="_gCSI" => "_ctrpv2")
common_gray_ctrpv2 = innerjoin(singleton_gray, singleton_ctrpv2, on=:pairs, renamecols="_gray" => "_ctrpv2")

### Common singletons for all three dataset
common_all = innerjoin(common_gCSI_ctrpv2, singleton_gray, on=:pairs)
rename!(common_all, :exp_id => :exp_id_gray)

### DO correlation
singletonCorrelation(common_gCSI_gray, ["gCSI", "gray"])
singletonCorrelation(common_gCSI_ctrpv2, ["gCSI", "ctrpv2"])
singletonCorrelation(common_gray_ctrpv2, ["gray", "ctrpv2"])