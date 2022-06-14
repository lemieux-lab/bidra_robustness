using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")

### Global Var
figure_prefix = "_generated_figures/discrepancies_replicates/"

function correlationBarPlot(df::DataFrame, t::String)
    p=Gadfly.plot(df, x=:description, y=:rₛ, color=:dataset, xgroup=:param, ygroup=:dataset,
                  Geom.subplot_grid( Geom.bar(position=:dodge)),
                  Guide.title(t))
    display(p)
end

function addMethods(df::DataFrame, method::String)
    df[!,:method] = repeat([method], nrow(df))
    return df
end

### Import correlation
ml_correlation = readCSV("_generated_data/mlCorrelations.csv", true, false, "")
median_correlation = readCSV("_generated_data/medianCorrelations.csv", true, false, "")
ci_correlation = readCSV("_generated_data/ciCorrelations.csv", true, false, "")
posterior_correlation = readCSV("_generated_data/posteriorCorrelations.csv", true, false, "")
qq_correlation = readCSV("_generated_data/qqCorrelations.csv", true, false, "")

### Plot correlations of individual methods
Gadfly.set_default_plot_size(5inch, 4inch)
correlationBarPlot(ml_correlation, "Marqaurdt-Levenberg Correlations")
correlationBarPlot(median_correlation, "BiDRA's median Correlations")
correlationBarPlot(qq_correlation, "BiDRA's QQ Correlations")

### Combine all correlations methods results and plot
ml_correlation = addMethods(ml_correlation, "ML")
median_correlation = addMethods(median_correlation, "Median")
qq_correlation = addMethods(qq_correlation, "QQ")
all_correlation = vcat(ml_correlation, median_correlation, qq_correlation)

Gadfly.set_default_plot_size(6inch, 4inch)
all_correlation_pairs = filter([:description, :param] => (x, y) -> x ∈ ["all pairs", "converged pairs"] && y != "aac", all_correlation)
p = Gadfly.plot(all_correlation_pairs, x=:method, y=:rₛ, color=:description, xgroup=:param, ygroup=:dataset, 
            Geom.subplot_grid(Geom.bar(position=:dodge)),
            Guide.title("All Pairs correlations"))
draw(PDF("$figure_prefix"*"all_dt_correlations.pdf", 6inch, 4inch), p)

### Correlations of random pairings
bidra_randomRep = readCSV("_generated_data/bidraRandomCorrelation.csv", true, false, "")

Gadfly.set_default_plot_size(10inch, 6inch)
p = Gadfly.plot(bidra_randomRep, x=:description, y=:rₛ, color=:dataset, xgroup=:param, ygroup=:method, 
            Geom.subplot_grid(Geom.boxplot()),
            Guide.title("Random Pairings Gray - R=$r"))
display(p)