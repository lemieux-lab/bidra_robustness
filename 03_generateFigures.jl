using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/utils.jl")

###### Global Var ####
figure_prefix = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_figures/03_qualitativeComparison/"

###### Data ##########
data_df = getRawData(["gCSI"], "/home/golem/scratch/labellec/_DATA/", true)
expId_list = unique(data_df.exp_id)

expIdSubset_list = ["NCI-H1648_AZ-628_8h", "Calu-1_PF-4708671_6b", "RERF-LC-MS_Gemcitabine_4b"]

gCSIml_df = getMLestimates(["gCSI"], false, missing)
gCSIml_df = filter(:exp_id => i -> i ∈ expId_list, gCSIml_df)
println(length(gCSIml_df.convergence[gCSIml_df.convergence]), " exp. that converge")

gCSIposterior_df = getBIDRAposterior("gCSI", expId_list)
gCSIposterior_df = innerjoin(gCSIposterior_df, gCSIml_df[:,[:exp_id, :convergence]], on=:exp_id)

println("Extreme values LDR: ML est. - BiDRA post.")
println("Max: ", maximum(gCSIml_df.LDR), "  ", maximum(gCSIposterior_df.LDR))
println("Minimum: ", minimum(gCSIml_df.LDR), "  ", minimum(gCSIposterior_df.LDR))

###### Grouping ##########
### Group ML estimate by convergence
filterLDR_number = filter(:LDR => x -> -50 ≤ x ≤ 150, gCSIml_df)
println(nrow(filterLDR_number), " exp. after LDR filtering based on numbers")
println(length(filterLDR_number.convergence[filterLDR_number.convergence]), " exp. that converge")
p = Gadfly.plot(filterLDR_number, x=:LDR, color=:convergence, Geom.histogram(bincount=200, density=true))

filterHDR_number = filter(:HDR => x -> -50 ≤ x ≤ 150, gCSIml_df)
println(nrow(filterHDR_number), " exp. after HDR filtering based on numbers")
println(length(filterHDR_number.convergence[filterHDR_number.convergence]), " exp. that converge")
p = Gadfly.plot(filterHDR_number, x=:HDR, color=:convergence, Geom.histogram(bincount=200, density=true))

filteIc50_number = filter(:ic50 => x -> -15 ≤ x ≤ 15, gCSIml_df)
println(nrow(filteIc50_number), " exp. after ic50 filtering based on numbers")
println(length(filteIc50_number.convergence[filteIc50_number.convergence]), " exp. that converge")

###### Median vs. ML estimations ##########
scaleColor = Scale.lab_gradient("gray95","black")

gCSImedian_df = combine(groupby(gCSIposterior_df, :exp_id), 
                        :HDR => median => :HDR_median,
                        :LDR => median => :LDR_median,
                        :ic50 => median => :ic50_median,
                        :exp_id => unique => :exp_id)
gCSImetrics_df = innerjoin(gCSImedian_df, gCSIml_df[:, [:LDR, :HDR, :ic50, :convergence, :exp_id]], on=:exp_id)
metricsSubset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, gCSImetrics_df)

## HDR comparison
λₖ = [0.4, 0.5, 0.1]
hdr_mm = MixtureModel([SkewNormal(0, 10, 1), Uniform(0, 100), SkewNormal(100, 20, -5)], λₖ)
hdr_prior = rand(hdr_mm, nrow(gCSIposterior_df))

metricsSubset = filter(:HDR => x -> -50 ≤ x ≤ 150, gCSImetrics_df)
correlationAnalysis(gCSImetrics_df.HDR_median, gCSImetrics_df.HDR)
correlationAnalysis(metricsSubset.HDR_median, metricsSubset.HDR)

Gadfly.set_default_plot_size(3inch, 2inch)
p1 = Gadfly.plot(metricsSubset, x=:HDR_median, y=:HDR, xintercept=[median(hdr_prior)],
                 Geom.vline(), Geom.hexbin(xbincount=80, ybincount=80),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.Cartesian(xmin=-20, ymin=-50, xmax=100, ymax=150),
                 Theme(panel_stroke="black"))
push!(p1, layer(metricsSubset_df, xintercept=:HDR_median, Geom.vline()))
push!(p1, layer(metricsSubset_df, yintercept=:HDR, Geom.hline()))
draw(PDF(figure_prefix*"med_mle_HDR.pdf", 3inch, 2inch), p1)

Gadfly.set_default_plot_size(3inch, 2inch)
p1x = Gadfly.plot(gCSIposterior_df, x=:HDR, color=:convergence, xintercept=[median(hdr_prior)],
                  Geom.vline(), Geom.histogram(density=true),
                  Coord.Cartesian(xmin=-20, xmax=100))
draw(PDF(figure_prefix*"med_mle_HDR_x_posterior.pdf", 3inch, 2inch), p1x)

Gadfly.set_default_plot_size(3inch, 2inch)
p1y = Gadfly.plot(metricsSubset, x=:HDR, color=:convergence,
                  Geom.histogram(density=true),
                  Coord.Cartesian(xmin=-50, xmax=150))
draw(PDF(figure_prefix*"med_mle_HDR_y_estimation.pdf", 3inch, 2inch), p1y)


## LDR comparison
ldr_prior = rand(Normal(100,10), 4000)

metricsSubset = filter(:LDR => x -> -50 ≤ x ≤ 150, gCSImetrics_df)
correlationAnalysis(gCSImetrics_df.LDR_median, gCSImetrics_df.LDR)
correlationAnalysis(metricsSubset.LDR_median, metricsSubset.LDR)

Gadfly.set_default_plot_size(3inch, 2inch)
p2 = Gadfly.plot(metricsSubset, x=:LDR_median, y=:LDR, xintercept=[median(ldr_prior)],
                 Geom.vline(), Geom.hexbin(xbincount=80, ybincount=80),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.Cartesian(xmin=65, ymin=50, xmax=150, ymax=150))
push!(p2, layer(metricsSubset_df, xintercept=:LDR_median, Geom.vline()))
push!(p2, layer(metricsSubset_df, yintercept=:LDR, Geom.hline()))
draw(PDF(figure_prefix*"med_mle_LDR.pdf", 3inch, 2inch), p2)

Gadfly.set_default_plot_size(3inch, 2inch)
p2x = Gadfly.plot(gCSIposterior_df, x=:LDR, color=:convergence, xintercept=[median(ldr_prior)],
                  Geom.vline(), Geom.histogram(density=true),
                  Coord.Cartesian(xmin=65, xmax=150))
draw(PDF(figure_prefix*"med_mle_LDR_x_posterior.pdf", 3inch, 2inch), p2x)

Gadfly.set_default_plot_size(3inch, 2inch)
p2y = Gadfly.plot(metricsSubset, x=:LDR, color=:convergence,
                  Geom.histogram(density=true),
                  Coord.Cartesian(xmin=50, xmax=150))
draw(PDF(figure_prefix*"med_mle_LDR_y_estimation.pdf", 3inch, 2inch), p2y)


## IC50 comparison
ic50_prior = rand(Normal(0,10), 4000)

metricsSubset = filter(:ic50 => x -> -15 ≤ x ≤ 15, gCSImetrics_df)
correlationAnalysis(gCSImetrics_df.ic50_median, gCSImetrics_df.ic50)
correlationAnalysis(metricsSubset.ic50_median, metricsSubset.ic50)

Gadfly.set_default_plot_size(3inch, 2inch)
p3 = Gadfly.plot(metricsSubset, x=:ic50_median, y=:ic50, xintercept=[median(ic50_prior)],
                 Geom.vline(), Geom.hexbin(xbincount=80, ybincount=80),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.Cartesian(ymin=-15, ymax=15))
push!(p3, layer(metricsSubset_df, xintercept=:ic50_median, Geom.vline()))
push!(p3, layer(metricsSubset_df, yintercept=:ic50, Geom.hline()))
draw(PDF(figure_prefix*"med_mle_ic50.pdf", 3inch, 2inch), p3)

Gadfly.set_default_plot_size(3inch, 2inch)
p3x = Gadfly.plot(gCSIposterior_df, x=:ic50, color=:convergence, xintercept=[median(ic50_prior)],
                  Geom.vline(), Geom.histogram(density=true),
                  Coord.Cartesian(xmin=-10, xmax=15))
draw(PDF(figure_prefix*"med_mle_ic50_x_posterior.pdf", 3inch, 2inch), p3x)

Gadfly.set_default_plot_size(3inch, 2inch)
p3y = Gadfly.plot(metricsSubset, x=:ic50, color=:convergence,
                  Geom.histogram(density=true),
                  Coord.Cartesian(xmin=-15, xmax=15))
draw(PDF(figure_prefix*"med_mle_ic5_y_estimation.pdf", 3inch, 2inch), p3y)

###### IC50 vs. experimental concentration ##########
### get standard deviation of raw data_df
concentrationBounds = percentile(data_df.Concentration, [2.5, 97.5])
gCSIdataMetrics_df = combine(groupby(data_df, :exp_id),
                            :Viability => std => :viability_std,
                            :Concentration => minimum => :concentration_min,
                            :Concentration => maximum => :concentration_max)
gCSImetrics_df = innerjoin(gCSImetrics_df, gCSIdataMetrics_df, on=:exp_id)

### compare ic50 estimation to std
Gadfly.set_default_plot_size(3inch, 2inch)
metricsSubset = filter(:ic50 => x -> -15 ≤ x ≤ 15, gCSImetrics_df)
p4 = Gadfly.plot(metricsSubset, x=:ic50, y=:viability_std,
                 Geom.hexbin(xbincount=80, ybincount=80),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.cartesian(xmin=-15, xmax=15, ymin=0),
                 Theme(panel_stroke= "black"))

expId_subset = filter(:exp_id => x -> x ∈ expIdSubset_list, metricsSubset)
push!(p4, layer(expId_subset, xintercept=:ic50, Geom.vline()))
push!(p4, layer(expId_subset, yintercept=:viability_std, Geom.hline()))
push!(p4, layer(xintercept=concentrationBounds, Geom.vline(color=["black"])))
draw(PDF(figure_prefix*"std_mle_ic50.pdf", 3inch, 2inch), p4)

### compare ic50 posterior to std
Gadfly.set_default_plot_size(3inch, 2inch)
gCSIposterior_df = innerjoin(gCSIposterior_df, gCSImetrics_df[:,[:exp_id, :viability_std]], on=:exp_id)
p5 = Gadfly.plot(gCSIposterior_df, x=:ic50, y=:viability_std,
                 Geom.hexbin(xbincount=100, ybincount=100),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.cartesian(xmin=-20, xmax=20, ymin=0),
                 Theme(panel_stroke="black"))

push!(p5, layer(expId_subset, yintercept=:viability_std, Geom.hline()))
push!(p5, layer(xintercept=concentrationBounds, Geom.vline(color=["black"])))

posteriorSubset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, gCSIposterior_df)
ic50_confidenceInterval_subset = combine(groupby(posteriorSubset_df, :exp_id), :ic50 => (x -> percentile(x, [2.5, 97.5])) => :ic50_ciBound)
push!(p5, layer(ic50_confidenceInterval_subset, xintercept=:ic50_ciBound, Geom.vline()))
draw(PDF(figure_prefix*"std_posterior_ic50.pdf", 3inch, 2inch), p5)

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