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


#### IC50 and HDR posterior
eff_metrics = [:ic50, :HDR]
bins_range = [collect(-5:0.1:5), collect(0:1:120)]

## prior Distributions
N = 100000
λₖ = [0.4, 0.5, 0.1]
mm = MixtureModel([SkewNormal(100, 10, 1), Uniform(0, 100), SkewNormal(0, 10, -1)], λₖ)
hdr_prior = rand(mm, N)
ic50_prior = rand(Normal(0, 10), N)

function plotPrior(dist::Array, plot_range::Vector)
    Gadfly.set_default_plot_size(3inch, 2inch)
    p1 = Gadfly.plot(x=dist, Geom.density(), Coord.cartesian(xmin=minimum(plot_range), xmax=maximum(plot_range)))
    p2 = Gadfly.plot(x=dist, Geom.histogram(), Coord.cartesian(xmin=minimum(plot_range), xmax=maximum(plot_range)))
    return p1, p2
end

p_ic50_prior = plotPrior(ic50_prior, bins_range[1]);
p_hdr_prior = plotPrior(hdr_prior, bins_range[2]);

fn = "ic50_prior_density.svg"
img = SVG(joinpath("figures", fn), 3inch, 2inch);
draw(img, p_ic50_prior[1]);

fn = "hdr_prior_density.svg"
img = SVG(joinpath("figures", fn), 3inch, 2inch);
draw(img, p_hdr_prior[1]);

## Create matrix
for p in 1:length(eff_metrics)
    pr = eff_metrics[p]
    pr_subset = posterior_data[:, [pr, :compound]]
    pr_posterior = reshape(pr_subset[:, pr], (4000, n))

    ## Calcute histogram binning
    pr_bins = bins_range[p]
    tmp_fit = mapreduce(r -> fit(Histogram, r, pr_bins).weights, vcat, eachcol(pr_posterior))
    y = repeat(compounds, inner=length(pr_bins)-1)
    x = repeat((pr_bins[1:end-1] + pr_bins[2:end]) ./ 2, outer=length(compounds))

    hist_long = DataFrame(x=x, y=y, count=tmp_fit)

    ## Plot
    scaleColor = Scale.lab_gradient("gray95","black")
    p = Gadfly.plot(hist_long, x=:x, y=:y, color=:count, 
                        Geom.rectbin(), 
                        Scale.color_continuous(colormap=scaleColor, minvalue=0),
                        Guide.title(string("$pr")),
                        Theme(panel_stroke= "black"));

    fn = string(string(pr), "_posterior.svg")
    img = SVG(joinpath("figures", fn), 4inch, 8inch);
    draw(img, p);
end


#### Compound subset: amongst the 20 first rank
prob = 0.9
topN = 20
tmp = counts_df_sorted[:, 1:topN]
prob_rank20 = DataFrame(prob_top=sum.(eachrow(tmp)), compound=compounds)

#### Select compound from top20
selected_prob = filter(:prob_top => x -> x >= prob, prob_rank20)
selected_posterior = filter(:compound => x -> x ∈ selected_prob.compound, posterior_data)

#### IC50 DAG
α = 0.05
selected_ic50 = selected_posterior[:, [:ic50, :compound]]
n = length(selected_prob.compound)

global dag = Graphs.DiGraph(n)

for i in 1:n
    for j in 1:n
        name_i = selected_prob.compound[i]
        name_j = selected_prob.compound[j]

        comp_i = filter(:compound => x -> x == name_i, selected_ic50)
        comp_j = filter(:compound => x -> x == name_j, selected_ic50)

        tmp = comp_i.ic50 .- comp_j.ic50
        Δ_prob = length(tmp[tmp .< 0]) / length(tmp)
        println(i, " -> ", j, ": ", Δ_prob)
        if Δ_prob >= (1 - α)
            Graphs.add_edge!(dag, i, j)
        end
    end
end

dag_reduce = Graphs.transitivereduction(dag)

p = graphplot(dag_reduce, names=selected_prob.compound, 
        method=:tree, fontsize=8, nodeshape=:rect, curves=false)
fn = "figures/ic50_dag.pdf"
Plots.pdf(p, fn)

### Compound selection probability
ldr_threshold = [-10, 10]
hdr_threshold = [50, 75] 

ldr_prob = combine(groupby(posterior_data, :compound), :LDR => sum∘(x -> (ldr_threshold[1] .< x .< ldr_threshold[2])) => :count_ldr)
ldr_prob[!, :prob_ldr] = ldr_prob.count_ldr ./ 4000
prob_select = innerjoin(prob_rank20, ldr_prob, on=:compound)

hdr_prob = combine(groupby(posterior_data, :compound), :HDR => sum∘(x -> hdr_threshold[1] .< x) => :count_hdr_50,
                                                    :HDR => sum∘(x -> hdr_threshold[2] .< x) => :count_hdr_80)
hdr_prob[!, :prob_hdr_50] = hdr_prob.count_hdr_50 ./ 4000
hdr_prob[!, :prob_hdr_80] = hdr_prob.count_hdr_80 ./ 4000
prob_select = innerjoin(prob_select, hdr_prob, on=:compound)

## sort based on ic50
prob_select = prob_select[indexin(compounds, prob_select.compound), :]
prob_select[!, :prob_select_50] = prob_select.prob_top .* prob_select.prob_ldr .* prob_select.prob_hdr_50
prob_select[!, :prob_select_80] = prob_select.prob_top .* prob_select.prob_ldr .* prob_select.prob_hdr_80
println(last(prob_select, 20))

## heatmap of prob
tmp = prob_select[:, [:prob_top, :prob_ldr, :prob_hdr_50, :prob_hdr_80, :prob_select_50, :prob_select_80, :compound]]
prob_long = DataFrames.stack(tmp, 1:6)

p = Gadfly.plot(prob_long, x=:variable, y=:compound, color=:value, 
                    Geom.rectbin(), 
                    Scale.color_continuous(colormap=scaleColor, minvalue=0,),
                    Guide.title(string("Selection_prob")),
                    Theme(panel_stroke= "black"));

fn = "prob_heatmap.svg"
img = SVG(joinpath("figures", fn), 4inch, 8inch);
draw(img, p);