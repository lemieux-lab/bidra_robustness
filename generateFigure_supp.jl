using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")

#### Correlation of experiments with multiple replicates
figure_prefix = "_generated_figures/supp_fig/multiRep_corr/"

## Sampling analysis
multiRep_ml = readCSV("_generated_data/multiRep_mlCorrelations.csv", true)
multiRep_ml[!, :method] = repeat(["Lev.-Mar."], nrow(multiRep_ml))

multiRep_median = readCSV("_generated_data/multiRep_medianCorrelations.csv", true)
multiRep_median[!, :method] = repeat(["Med. Posterior"], nrow(multiRep_median))

multiRep_qq = readCSV("_generated_data/multiRep_qqCorrelations.csv", true)
multiRep_qq[!, :method] = repeat(["QQ Posterior"], nrow(multiRep_qq))

multiRep_all = vcat(multiRep_ml, multiRep_median, multiRep_qq)

function boxplotCorr(df::DataFrame, T::String, coef::Symbol)
    sg = Geom.subplot_grid(Geom.boxplot())
    sg.coord = Coord.cartesian(ymin=0, ymax=1.)
    p = Gadfly.plot(df, x=:param, y=coef, color=:param, xgroup=:dataset, 
                    sg,
                    Guide.title(T))
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