using DataFrames, CSV
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots, Plots
#using Cairo, Fontconfig
using Graphs, GraphRecipes, LightGraphs

include("utils_iric.jl")

## Get data
posterior_path = "results/"
posterior_data = getPosterior(posterior_path)
viability_data = getViability()
compounds = unique(posterior_data.compound)
n = length(compounds)

## Get matrice
eff_metrics = :ic50             ## main metrics
rank_order = false              ## rev. order

pr_subset = posterior_data[:, [eff_metrics, :compound]]
median_data = combine(groupby(posterior_data, :compound), eff_metrics => median => :median)
sort!(median_data, :median, rev=!rank_order)
## row : MCMC sampling
## col : compound
values_matrix = reshape(pr_subset[:, eff_metrics], (4000, n))
ranks_matrix = mapslices(x -> ordinalrank(x, rev=rank_order), values_matrix, dims=2)
counts_rank = countmap.(eachcol(ranks_matrix))
## row : compounds
## col : prob
counts_matrix = zeros(Float64, n, n)
for i in 1:n
    tmp = counts_rank[i]

    for (index, value) in pairs(tmp)
        counts_matrix[i,index] = value/4000
    end
end

## Define ranks-prob df
ranks_label = "R" .* numstring.([i for i in 1:n], 4)
counts_df = DataFrame(counts_matrix, :auto)
rename!(counts_df, ranks_label)
counts_df[!, :compound] = compounds

## order based on IC50 median values
counts_df_sorted = counts_df[indexin(median_data.compound, counts_df.compound), :]
counts_df_long = DataFrames.stack(counts_df_sorted, 1:n)

## update ordering of posterior and compounds
compounds = counts_df_sorted.compound
posterior_data[!, :sortby] = indexin(posterior_data.compound, compounds)
sort!(posterior_data, :sortby)

scaleColor = Scale.lab_gradient("gray95","black")
Gadfly.set_default_plot_size(10inch, 8inch)
p = Gadfly.plot(counts_df_long, x=:variable, y=:compound, color=:value, 
                Geom.rectbin(), 
                Scale.color_continuous(colormap=scaleColor, minvalue=0, maxvalue=1),
                #Guide.title(string("$sample_name - $eff_metrics")),
                );
fn = string(eff_metrics, "_ranking.svg")
img = SVG(joinpath("figures", fn), 8inch, 8inch);
draw(img, p);


#### SD of viability responses
sd_df = combine(groupby(viability_data, [:anonym_id]), :inhib_col => std => :viab_std)

p = Gadfly.plot(sd_df[indexin(compounds, sd_df.anonym_id), :], 
                x=:anonym_id, y=:viab_std, yintercept=[20],
                Geom.bar(), Geom.hline());

fn = string(eff_metrics, "_SD_viab.svg")
img = SVG(joinpath("figures", fn), 8inch, 3inch);
draw(img, p);
