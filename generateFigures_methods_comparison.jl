using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")
include("plot_utils.jl")

###### Global Var ####
figure_prefix = "_generated_figures/methods_comparison/"
datasets = ["gCSI", "gray", "ctrpv2"]
eff_metrics = [:HDR, :LDR, :ic50, :slope, :aac]
metrics_bounds = Dict(:HDR=>[-50,150], :LDR=>[70,150], :ic50=>[-10,10], :slope=>[0,10], :aac=>[0,100])
expIdSubset_list = ["NCI-H1648_AZ-628_8h", "Calu-1_PF-4708671_6b", "RERF-LC-MS_Gemcitabine_4b", "HCC78_Lapatinib_11a"]
scaleColor = Scale.lab_gradient("gray95","black")
N = 1000

getPriors = Dict{String, Function}()
getPriors["HDR"] = (N) -> get_HDR_prior(N)
getPriors["LDR"] = (N) -> get_LDR_prior(N)
getPriors["ic50"] = (N) -> get_ic50_prior(N)
getPriors["slope"] = (N) -> get_slope_prior(N)
getPriors["aac"] = (N) -> [0,0]

###### Data ##########
i = 2

dt = datasets[2]
println("**********************", dt, "**********************")
println("1. Get data and metrics")
data_df = getRawData([dt], "/home/golem/scratch/labellec/_DATA/", true)
expId_list = unique(data_df.exp_id)

ml_df = getMLestimates([dt], false, missing)
tmp = get_converged(ml_df)

posterior_df = getBIDRAposterior(dt, expId_list)
#posterior_df = innerjoin(posterior_df, ml_df[:,[:exp_id, :convergence]], on=:exp_id)


###### Median vs. ML estimations ##########
println("2. Comparing LM estimate to posterior median")
metrics_df = get_median(posterior_df, eff_metrics, ml_df)

for em in eff_metrics
    println("---> ", em)
    println("------> Plotting hexbin")
    em_subset_df = metricZoom_subset(metrics_df, em, metrics_bounds[em][1], metrics_bounds[em][2])
    em_prior = getPriors[string(em)](N)
    p = medianML_hexbin_plot(em_subset_df, em, metrics_bounds[em][1], metrics_bounds[em][2], em_prior)

    if dt == datasets[1]
        exp_subset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, metrics_df)
        push!(p, layer(exp_subset_df, x=Symbol(String(em)*"_median"), y=em, Geom.point()))
    end

    draw(PDF(figure_prefix*dt*"_med_mle_"*string(em)*".pdf", 3inch, 2inch), p)

    println("------> Plotting contour")
    em_mtx = kde((em_subset_df[:, Symbol(String(em)*"_median")], em_subset_df[:, em]))
    p = medianML_contour_plot(em_mtx, metrics_bounds[em][1], metrics_bounds[em][2], em_prior)
    draw(PDF(figure_prefix*dt*"_med_mle_"*string(em)*"_contour.pdf", 3inch, 2inch), p)
end


###### IC50 vs. experimental concentration ##########
println("3. Looking @ IC50 and experimental dose range")
concentrationBounds, metrics_df = get_data_metrics(data_df, metrics_df)
posterior_df = innerjoin(posterior_df, metrics_df[:,[:exp_id, :viability_std]], on=:exp_id)

### compare ic50 estimation to std
println("---> Plotting IC50 est. vs. std")
subset_df = metricZoom_subset(metrics_df, :ic50, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2])
p = ic50_std_plot(subset_df, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2], concentrationBounds)

if dt == datasets[1]
    exp_subset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, subset_df)
    push!(p, layer(exp_subset_df, x=:ic50, y=:viability_std, Geom.point()))
end

draw(PDF(figure_prefix*dt*"_std_mle_ic50.pdf", 3inch, 2inch), p)

println("---> Plotting IC50 est. vs. std contour")
mtx = kde((subset_df[:, :ic50], subset_df[:, :viability_std]))
p = ic50_std_contour_plot(mtx, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2], concentrationBounds)
draw(PDF(figure_prefix*dt*"_std_mle_ic50_contour.pdf", 3inch, 2inch), p)

println("---> Plotting IC50 posterior vs. std")
subset_df = metricZoom_subset(posterior_df, :ic50, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2])
p = ic50_std_plot(subset_df, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2], concentrationBounds)
draw(PDF(figure_prefix*dt*"_std_posterior_ic50.pdf", 3inch, 2inch), p)

println("---> Plotting IC50 posterior vs. std contour")
mtx = kde((subset_df[:, :ic50], subset_df[:, :viability_std]))
p = ic50_std_contour_plot(mtx, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2], concentrationBounds)
draw(PDF(figure_prefix*dt*"_std_posterior_ic50_contour.pdf", 3inch, 2inch), p)

println("---> Calculating Prob. of IC50 being outside of experimental dose range")
withinDose_prob, posterior_df
### compare std and P(max < ic50)
gCSIposterior_df = innerjoin(gCSIposterior_df, gCSImetrics_df[:,[:exp_id, :concentration_min, :concentration_max]], on=:exp_id)
gCSIposterior_df[:, :withinDose] = gCSIposterior_df.concentration_min .< gCSIposterior_df.ic50 .< gCSIposterior_df.concentration_max

gCSImetrics_df[:, :withinDoseML] = gCSImetrics_df.concentration_min .< gCSImetrics_df.ic50 .< gCSImetrics_df.concentration_max
gCSImetrics_df[:, :outsideDoseML] = .!gCSImetrics_df.withinDoseML

withinDose_prob = combine(groupby(gCSIposterior_df, :exp_id), :withinDose => sum => :withinCount,
                                                              :viability_std => unique => :viability_std)
withinDose_prob[:, :withinProb] = withinDose_prob.withinCount ./ 4000
withinDose_prob[:, :outsideProb] = 1. .- withinDose_prob.withinProb
withinDose_prob = innerjoin(withinDose_prob, gCSImetrics_df[:, [:exp_id, :convergence, :outsideDoseML]], on=:exp_id)

Gadfly.set_default_plot_size(5inch, 2inch)
p6 = Gadfly.plot(withinDose_prob, x=:outsideProb, y=:viability_std, 
                 Geom.hexbin(xbincount=120, ybincount=90),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.cartesian(xmin=0, ymin=0, xmax=1),
                 Theme(panel_stroke="black"))

correlationAnalysis(withinDose_prob.outsideProb, withinDose_prob.viability_std)
withinDoseSubset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, withinDose_prob)
push!(p6, layer(withinDoseSubset_df, xintercept=:outsideProb, Geom.vline()))
push!(p6, layer(withinDoseSubset_df, yintercept=:viability_std, Geom.hline()))
draw(PDF(figure_prefix*"std_prob_ic50.pdf", 5inch, 2inch), p6)

Gadfly.set_default_plot_size(2inch, 2inch)
p6xa = Gadfly.plot(withinDose_prob, x=:outsideDoseML, y=:outsideProb, 
                   Geom.boxplot(),)
draw(PDF(figure_prefix*"boxplot_ml_prob_ic50.pdf", 2inch, 2inch), p6xa)
p6xb = Gadfly.plot(withinDose_prob, x=:convergence, y=:outsideProb, 
                   Geom.boxplot(),)
draw(PDF(figure_prefix*"boxplot_conv_prob_ic50.pdf", 2inch, 2inch), p6xb)

Gadfly.set_default_plot_size(3inch, 2inch)
p6ya = Gadfly.plot(withinDose_prob, x=:viability_std, color=:convergence, 
                   Geom.histogram(density=true),
                   Coord.cartesian(xmin=0, xmax=60))
draw(PDF(figure_prefix*"hist_conv_viabStd_ic50.pdf", 3inch, 2inch), p6ya)

p6yb = Gadfly.plot(withinDose_prob, x=:viability_std, color=:outsideDoseML,
                   Geom.histogram(density=true),
                   Coord.cartesian(xmin=0, xmax=60))
draw(PDF(figure_prefix*"hist_ml_viabStd_ic50.pdf", 3inch, 2inch), p6yb)


###### MAD vs. ML estimations ##########
gCSImad_df = combine(groupby(gCSIposterior_df, :exp_id), 
                        :HDR => mad => :HDR_mad,
                        :LDR => mad => :LDR_mad,
                        :ic50 => mad => :ic50_mad)
filter(:exp_id => x -> x ∈ expIdSubset_list, gCSImad_df)

gCSImetrics_df = innerjoin(gCSImetrics_df, gCSImad_df, on=:exp_id)

## HDR comparison
Gadfly.set_default_plot_size(3inch, 2inch)
metricsSubset = filter(:HDR => x -> -50 ≤ x ≤ 150, gCSImetrics_df)
p7 = Gadfly.plot(metricsSubset, x=:HDR, y=:HDR_mad, 
                 Geom.hexbin(xbincount=120, ybincount=100),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.cartesian(xmin=-50, xmax=150, ymin=0),
                 Theme(panel_stroke="black"))

expIDsubset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, gCSImetrics_df)
push!(p7, layer(expIDsubset_df, xintercept=:HDR, Geom.vline()))
push!(p7, layer(expIDsubset_df, yintercept=:HDR_mad, Geom.hline()))
draw(PDF(figure_prefix*"mad_hdr.pdf", 3inch, 2inch), p7)

###### postDif vs. ML estimations ##########
gCSIposteriorDiff_df = combine(groupby(gCSIposterior_df, :exp_id), 
                                        :HDR => (x -> percentile(x, 97.5) - percentile(x, 2.5)) => :HDR_postDiff,
                                        :LDR => (x -> percentile(x, 97.5) - percentile(x, 2.5)) => :LDR_postDiff,
                                        :ic50 => (x -> percentile(x, 97.5) - percentile(x, 2.5)) => :ic50_postDiff)
filter(:exp_id => x -> x ∈ expIdSubset_list, gCSIposteriorDiff_df)

gCSImetrics_df = innerjoin(gCSImetrics_df, gCSIposteriorDiff_df, on=:exp_id)

## HDR comparison
Gadfly.set_default_plot_size(3inch, 2inch)
metricsSubset = filter(:HDR => x -> -50 ≤ x ≤ 150, gCSImetrics_df)
p8 = Gadfly.plot(metricsSubset, x=:HDR, y=:HDR_postDiff, 
                 Geom.hexbin(xbincount=120, ybincount=100),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.cartesian(xmin=-50, xmax=150, ymin=0),
                 Theme(panel_stroke="black"))

expIDsubset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, gCSImetrics_df)
push!(p8, layer(expIDsubset_df, xintercept=:HDR, Geom.vline()))
push!(p8, layer(expIDsubset_df, yintercept=:HDR_postDiff, Geom.hline()))
draw(PDF(figure_prefix*"postDiff_hdr.pdf", 3inch, 2inch), p8)

## IC50 comparison
Gadfly.set_default_plot_size(3inch, 2inch)
metricsSubset = filter(:ic50 => x -> -15 ≤ x ≤ 15, gCSImetrics_df)
p9 = Gadfly.plot(metricsSubset, x=:ic50, y=:ic50_postDiff, 
                 Geom.hexbin(xbincount=120, ybincount=100),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.cartesian(xmin=-15, xmax=15, ymin=0),
                 Theme(panel_stroke="black"))

expIDsubset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, gCSImetrics_df)
push!(p9, layer(expIDsubset_df, xintercept=:ic50, Geom.vline()))
push!(p9, layer(expIDsubset_df, yintercept=:ic50_postDiff, Geom.hline()))
draw(PDF(figure_prefix*"postDiff_ic50.pdf", 3inch, 2inch), p9)

#### Example #4
## exp_id = HCC78_Lapatinib_11a, KP-3_Erastin_4b, KMM-1_Dabrafenib_11d
e = "HCC78_Lapatinib_11a"
Gadfly.set_default_plot_size(4inch, 2inch)
viab_subset = filter(:exp_id => x -> x == e, data_df)
p10 = Gadfly.plot(layer(viab_subset, x=:Concentration, y=:Viability, color=:exp_id, Geom.point()),
                  Coord.cartesian(ymin=-10, ymax=110),
                  Theme(panel_stroke="black"))

subset_posterior = getBIDRAposterior("gCSI", [e])
postCurves = getPosteriorCurves(subset_posterior, -4, 1)
xDose = parse.(Float64, names(postCurves))

med_curve = median.(eachcol(postCurves))
upper_curve = percentile.(eachcol(postCurves), 97.5)
lower_curve = percentile.(eachcol(postCurves), 2.5)
push!(p10, layer(x=xDose, y=med_curve, color=[e], linestyle=[:dot], Geom.line()))
        
mleFit = llogistic(filter(:exp_id => x -> x == e, gCSIml_df)[1, [:LDR, :HDR, :ic50, :slope]])
yViability = mleFit.(xDose)
push!(p10, layer(x=xDose, y=yViability, color=[e], Geom.line()))

push!(p10, layer(x=xDose, ymin=lower_curve, ymax=upper_curve, color=[e], alpha=[0.3], Geom.ribbon()))
draw(PDF(figure_prefix*"curveExample.pdf", 4inch, 2inch), p10)

Gadfly.set_default_plot_size(3inch, 2inch)
x = percentile(subset_posterior.HDR, [2.5, 97.5])
ymin = [0.,0.]
ymax = ymin .+ 0.04
hdr = filter(:exp_id => x -> x == e, gCSIml_df)[:, :HDR]
p11 = Gadfly.plot(subset_posterior, x=:HDR, Geom.histogram(density=true),
                  layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon()),
                  layer(xintercept=hdr, Geom.vline()),
                  Theme(panel_stroke="black"))
draw(PDF(figure_prefix*"posteriorExample_HDR.pdf", 3inch, 2inch), p11)

Gadfly.set_default_plot_size(3inch, 2inch)
x = percentile(subset_posterior.ic50, [2.5, 97.5])
ymin = [0.,0.]
ymax = ymin .+ 4
ic50 = filter(:exp_id => x -> x == e, gCSIml_df)[:, :ic50]
p12 = Gadfly.plot(subset_posterior, x=:ic50, Geom.histogram(density=true),
                  layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon()),
                  layer(xintercept=ic50, Geom.vline()),
                  Theme(panel_stroke="black"))
draw(PDF(figure_prefix*"posteriorExample_IC50.pdf", 3inch, 2inch), p12)