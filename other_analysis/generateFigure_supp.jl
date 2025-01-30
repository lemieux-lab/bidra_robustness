using HypothesisTests
using Gadfly
using Cairo, Fontconfig
using DataFrames, HDF5, JLD2

include("../utils.jl")
include("../plot_utils.jl")

#### Correlation of experiments with multiple replicates
figure_prefix = "_generated_figures/supp_fig/multiRep_corr"

## Numbers of exp./pairs
data_prefix = "public_datasets/"
for dt in ["gray"]#, "gCSI", "ctrpv2"]
    grouped_expID = load(joinpath(data_prefix, "rep_more2_pairing.jld2"))[dt]
    expID_by_pairs = [length(g.pairs) for g in  groupby(grouped_expID, :pairs)]

    tmp = countmap(expID_by_pairs)
    p = Gadfly.plot(x = collect(keys(tmp)), y= collect(values(tmp)), Geom.bar(),)
    draw(PDF(figure_prefix*"_$dt"*"_count_exp.pdf", 6inch, 4inch), p)
end

## Sampling analysis
multiRep_ml = readCSV("_generated_data/multiRep_mlCorrelations.csv", true);
multiRep_ml[!, :method] = repeat(["Lev.-Mar."], nrow(multiRep_ml))

multiRep_median = readCSV("_generated_data/multiRep_medianCorrelations.csv", true);
multiRep_median[!, :method] = repeat(["Med. Posterior"], nrow(multiRep_median))

multiRep_qq = readCSV("_generated_data/multiRep_qqCorrelations.csv", true);
multiRep_qq[!, :method] = repeat(["QQ Posterior"], nrow(multiRep_qq));

multiRep_all = vcat(multiRep_ml, multiRep_median, multiRep_qq);

function boxplotCorr(df::DataFrame, T::String, coef::Symbol)
    sg = Geom.subplot_grid(Geom.boxplot())
    sg.coord = Coord.cartesian(ymin=0, ymax=1.)
    p = Gadfly.plot(df, x=:param, y=coef, color=:param, xgroup=:dataset, 
                    sg,
                    Guide.title(T));
    return p
end

p = boxplotCorr(multiRep_ml, "10,000 Rep. ML Correlation", :rₛ)
draw(PDF(figure_prefix*"sampling_rep_ml_corr_spearman.pdf", 10inch, 6inch), p)

p = boxplotCorr(multiRep_ml, "10,000 Rep. ML Correlation", :r)
draw(PDF(figure_prefix*"sampling_rep_ml_corr_pearson.pdf", 10inch, 6inch), p)

p = boxplotCorr(multiRep_median, "10,000 Rep. median Correlation", :rₛ)
draw(PDF(figure_prefix*"sampling_rep_median_corr_spearmann.pdf", 10inch, 6inch), p)

p = boxplotCorr(multiRep_median, "10,000 Rep. median Correlation", :r)
draw(PDF(figure_prefix*"sampling_rep_median_corr_pearson.pdf", 10inch, 6inch), p)

p = boxplotCorr(multiRep_qq, "10,000 Rep. QQ Correlation", :rₛ)
draw(PDF(figure_prefix*"sampling_rep_qq_corr_spearman.pdf", 10inch, 6inch), p)

p = boxplotCorr(multiRep_qq, "10,000 Rep. QQ Correlation", :r)
draw(PDF(figure_prefix*"sampling_rep_qq_corr_pearson.pdf", 10inch, 6inch), p)

p = Gadfly.plot(multiRep_all, x=:method, y=:rₛ, color=:method, xgroup=:param,ygroup=:dataset,
                Geom.subplot_grid(Geom.boxplot),
                Guide.title("Correlation of 10,000 Rep. Sampling for Multi. Rep. expeirments"))
draw(PDF(figure_prefix*"sampling_rep_all_spearman.pdf", 10inch, 10inch), p)

p = Gadfly.plot(multiRep_all, x=:method, y=:r, color=:method, xgroup=:param,ygroup=:dataset,
                Geom.subplot_grid(Geom.boxplot),
                Guide.title("Correlation of 10,000 Rep. Sampling for Multi. Rep. expeirments"))
draw(PDF(figure_prefix*"sampling_rep_all_pearson.pdf", 10inch, 10inch), p)


#### Correlation of experiments with multiple replicates
figure_prefix_across = "_generated_figures/supp_fig/across_dataset/"

## Correlation analysis
acrossCorr = readCSV("_generated_data/acrossDataset_correlations.csv", true)
p = Gadfly.plot(acrossCorr, x=:method, y=:rₛ, xgroup=:param, ygroup=:pair, color=:method,
                Geom.subplot_grid(Geom.bar()))
draw(PDF(figure_prefix_across*"correlation_by_method_spearman.pdf", 6inch, 6inch), p)

p = Gadfly.plot(acrossCorr, x=:method, y=:r, xgroup=:param, ygroup=:pair, color=:method,
                Geom.subplot_grid(Geom.bar()))
draw(PDF(figure_prefix_across*"correlation_by_method_pearson.pdf", 6inch, 6inch), p)

#### How to assess if we can trust the results
figure_prefix_trust = "_generated_figures/supp_fig/trust_metrics/"
metrics_bounds = Dict(:HDR=>[-50,150], :LDR=>[70,150], :ic50=>[-10,10], :slope=>[0,10]);


for dt in ["gray", "gCSI", "ctrpv2"]
    exp_id = getExpId_h5(dt)

    data_df = getRawData_h5(dt, false)

    ## Response at max concentration
    tmp = combine(groupby(data_df, [:exp_id, :Concentration]), :Viability => mean => :avg_viability)
    metrics_df = combine(groupby(tmp, :exp_id), [:Concentration, :avg_viability] => ((a,b) -> b[argmax(a)]) => :avg_viab_max_conc)

    ## Response SD
    tmp = combine(groupby(data_df, :exp_id), :Viability => std => :sd_viab)
    metrics_df = innerjoin(metrics_df, tmp, on=:exp_id)

    ## LM convergence status + Metrics
    tmp = getMLestimates(dt)
    metrics_df = innerjoin(metrics_df, tmp[:, [:exp_id, :rmsd, :convergence, :LDR, :HDR, :ic50, :slope]], on=:exp_id)

    ## ΔHDR, ΔIC50, ΔSlope
    posterior = getPosterior_h5(dt, false, exp_id)
    metrics_df = get_posterior_diff(posterior, metrics_df)

    #### LM metrics
    groupX = filter(:convergence => x -> x == true, metrics_df)
    groupY = filter(:convergence => x -> x == false, metrics_df)

    # Viability of max concentration
    p = Gadfly.plot(metrics_df, x=:avg_viab_max_conc, color=:convergence, Geom.histogram());
    draw(PDF(figure_prefix_trust*"$dt"*"_distributions_avg_viab_max_conc.pdf", 4inch, 4inch), p);
    utest = pvalue(MannWhitneyUTest(Float64.(groupX.avg_viab_max_conc), Float64.(groupY.avg_viab_max_conc)), tail=:both)
    p = Gadfly.plot(metrics_df, x=:convergence, y=:avg_viab_max_conc, Geom.boxplot(),
                    Guide.title("MWU-test: $utest"));
    draw(PDF(figure_prefix_trust*"$dt"*"_boxplot_avg_viab_max_conc.pdf", 4inch, 4inch), p)      

    # SD of viability
    p = Gadfly.plot(metrics_df, x=:sd_viab, color=:convergence, Geom.histogram());
    draw(PDF(figure_prefix_trust*"$dt"*"_distributions_sd_viab.pdf", 4inch, 4inch), p);
    utest = pvalue(MannWhitneyUTest(Float64.(groupX.sd_viab), Float64.(groupY.sd_viab)), tail=:both)
    p = Gadfly.plot(metrics_df, x=:convergence, y=:sd_viab, yintercept=[20], Geom.boxplot(), Geom.hline(),
                    Guide.title("MWU-test: $utest"));
    draw(PDF(figure_prefix_trust*"$dt"*"_boxplot_sd_viab.pdf", 4inch, 4inch), p) 

    # RMSE of LM
    p = Gadfly.plot(metrics_df, x=:rmsd, color=:convergence, Geom.histogram());
    draw(PDF(figure_prefix_trust*"$dt"*"_distributions_rmsd.pdf", 4inch, 4inch), p);
    utest = pvalue(MannWhitneyUTest(Float64.(groupX.rmsd), Float64.(groupY.rmsd)), tail=:both)
    p = Gadfly.plot(metrics_df, x=:convergence, y=:rmsd, Geom.boxplot(),
                    Guide.title("MWU-test: $utest"));
    draw(PDF(figure_prefix_trust*"$dt"*"_boxplot_rmsd.pdf", 4inch, 4inch), p)

    # Metrics vs. RMSE
    tmp_rmsd = metrics_df[:, [:exp_id, :rmsd, :convergence, :LDR, :HDR, :ic50, :slope]]
    tmp_stack = stack(tmp_rmsd, [:LDR, :HDR, :ic50, :slope])
    tmp_filter = DataFrame()

    for pr in keys(metrics_bounds)
        tmp = filter([:variable, :value] => (a, b) -> (a == String(pr) && metrics_bounds[pr][1] ≤ b ≤ metrics_bounds[pr][2]), tmp_stack)
        tmp_filter = vcat(tmp_filter, tmp)
    end
    p = Gadfly.plot(tmp_filter, x=:value, y=:rmsd, xgroup=:variable, color=:convergence, Geom.subplot_grid(Geom.point(), free_x_axis=true));
    draw(PDF(figure_prefix_trust*"$dt"*"_metrics_rmsd.pdf", 10inch, 4inch), p)
    println("metrics vs. RMSD")
    print(countmap(tmp_filter.variable))
    println()

    # Metrics vs. SD_viab
    tmp_sd = metrics_df[:, [:exp_id, :sd_viab, :convergence, :LDR, :HDR, :ic50, :slope]]
    tmp_stack = stack(tmp_sd, [:LDR, :HDR, :ic50, :slope])
    tmp_filter = DataFrame()

    for pr in keys(metrics_bounds)
        tmp = filter([:variable, :value] => (a, b) -> (a == String(pr) && metrics_bounds[pr][1] ≤ b ≤ metrics_bounds[pr][2]), tmp_stack)
        tmp_filter = vcat(tmp_filter, tmp)
    end

    scaleColor = Scale.lab_gradient("gray95","black")
    println(combine(groupby(tmp_filter, :variable), nrow))
    p = Gadfly.plot(tmp_filter, x=:value, y=:sd_viab, xgroup=:variable, yintercept=[20],
                    Geom.subplot_grid(Geom.hexbin(xbincount=80, ybincount=80), Geom.hline(), free_x_axis=true),
                    Scale.color_continuous(colormap=scaleColor, minvalue=1));
    draw(PDF(figure_prefix_trust*"$dt"*"_metrics_sd_viab_hexbin.pdf", 10inch, 4inch), p)

    p = Gadfly.plot(tmp_filter, x=:value, y=:sd_viab, xgroup=:variable, color=:convergence, yintercept=[20],
                    Geom.subplot_grid(Geom.point(), Geom.hline(), free_x_axis=true));
    draw(PDF(figure_prefix_trust*"$dt"*"_metrics_sd_viab.pdf", 10inch, 4inch), p)
    println("metrics vs. SD viab")
    print(countmap(tmp_filter.variable))
    println()


    #### Posterior metrics
    ## Comparing delta to SD-viab
   
    tmp_sd = metrics_df[:, [:exp_id, :sd_viab, :convergence, :LDR_postDiff, :HDR_postDiff, :ic50_postDiff, :slope_postDiff]]
    tmp_stack = stack(tmp_sd, [:LDR_postDiff, :HDR_postDiff, :ic50_postDiff, :slope_postDiff])
    p = Gadfly.plot(tmp_stack, x=:value, y=:sd_viab, yintercept=[20], xgroup=:variable,
                    Geom.subplot_grid(Geom.hexbin(xbincount=80, ybincount=80), Geom.hline, free_x_axis=true),
                    Scale.color_continuous(colormap=scaleColor, minvalue=1));
    draw(PDF(figure_prefix_trust*"$dt"*"_Δmetrics_sd_viab.pdf", 10inch, 4inch), p)
    println("Δmetrics vs. SD viab")
    print(countmap(tmp_stack.variable))
    println()


    tmp_posterior = innerjoin(posterior[:, [:exp_id, :LDR, :HDR, :ic50, :slope]], metrics_df[:, [:exp_id, :sd_viab]], on=:exp_id)
    tmp_stack = stack(tmp_posterior, [:LDR, :HDR, :ic50, :slope])
    p = Gadfly.plot(tmp_stack, x=:value, y=:sd_viab, yintercept=[20], xgroup=:variable,
                    Geom.subplot_grid(Geom.hexbin(xbincount=80, ybincount=80), Geom.hline, free_x_axis=true),
                    Scale.color_continuous(colormap=scaleColor, minvalue=1));
    draw(PDF(figure_prefix_trust*"$dt"*"_posterior_metrics_sd_viab.pdf", 10inch, 4inch), p)


    #p = Gadfly.plot(tmp_stack, x=:value, y=:sd_viab, yintercept=[20], color=:convergence, xgroup=:variable,
    #                Geom.subplot_grid(Geom.point(), Geom.hline, free_x_axis=true));
    #draw(PDF(figure_prefix_trust*"$dt"*"_Δmetrics_sd_viab_convergence.pdf", 10inch, 4inch), p)
end
