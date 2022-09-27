using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")

#### Correlation of expeirments with multiple replicate
figure_prefix = "_generated_figures/supp_fig/multiRep_corr/"

## Sampling analysis
multiRep_ml = readCSV("_generated_data/multiRep_mlCorrelations.csv", true)
multiRep_median = readCSV("_generated_data/multiRep_medianCorrelations.csv", true)
multiRep_qq = readCSV("_generated_data/multiRep_qqCorrelations.csv", true)

function boxplotCorr(df::DataFrame, T::String)
    sg = Geom.subplot_grid(Geom.boxplot())
    sg.coord = Coord.cartesian(ymin=0, ymax=1.)
    p = Gadfly.plot(df, x=:dataset, y=:rₛ, color=:dataset, xgroup=:param, 
                    sg,
                    Guide.title(T))
    return p
end

p = boxplotCorr(multiRep_ml, "10,000 Rep. ML Correlation")
draw(PDF(figure_prefix*"sampling_rep_ml_corr.pdf", 10inch, 6inch), p)

p = boxplotCorr(multiRep_median, "10,000 Rep. median Correlation")
draw(PDF(figure_prefix*"sampling_rep_median_corr.pdf", 10inch, 6inch), p)

p = boxplotCorr(multiRep_qq, "10,000 Rep. QQ Correlation")
draw(PDF(figure_prefix*"sampling_rep_qq_corr.pdf", 10inch, 6inch), p)



## Metrics SD by grouped of multi replicates
dt = "gray"