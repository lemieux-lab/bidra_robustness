using DataFrames, DataStructures
using Distributions, Statistics, StatsBase, KernelDensity
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")

### Given posterior distribution, calculate the median values
### for all efficacy metrics, and create a metric df 
function get_median(posterior_df, eff_metrics, ml_df)
    median_df = combine(groupby(posterior_df, :exp_id), 
                    eff_metrics[1] => median,
                    eff_metrics[2] => median,
                    eff_metrics[3] => median,
                    eff_metrics[4] => median,
                    eff_metrics[5] => median,
                    :exp_id => unique => :exp_id)
    metrics_df = innerjoin(median_df, ml_df[:, vcat(eff_metrics, [:convergence, :exp_id])], on=:exp_id)
    return metrics_df
end

### Select experiments whose metrics values are within a 
### given range
function metricZoom_subset(metrics_df, mt, lb, ub)
    tmp = filter(mt => x -> lb ≤ x ≤ ub, metrics_df)
    println("------> ", length(unique(tmp.exp_id)), " experiments out of ", length(unique(metrics_df.exp_id)), " for ", mt)
    tmp[!,mt] = convert.(Float64, tmp[!,mt])
    return tmp
end


### Plot estimates vs. median with hexbin for a given range of of values
function medianML_hexbin_plot(df, mt, lb, ub, prior)
    Gadfly.set_default_plot_size(3inch, 2inch)
    p = Gadfly.plot(df, x=Symbol(String(mt)*"_median"), y=mt, xintercept=[median(prior)],
                 Geom.vline(), Geom.hexbin(xbincount=80, ybincount=80),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.Cartesian(xmin=lb, ymin=lb, xmax=ub, ymax=ub),
                 Theme(panel_stroke="black"),
                 Guide.title("N="*string(nrow(df))))
    return p
end

### Plot estimates vs. median with contour for a given range of of values
function medianML_contour_plot(mtx, lb, ub, prior)
    Gadfly.set_default_plot_size(3inch, 2inch)
    p = Gadfly.plot(z=mtx.density, x=mtx.x, y=mtx.y, xintercept=[median(prior)],
                 Geom.vline(), Geom.contour(levels=8),
                 Scale.color_continuous(colormap=scaleColor, minvalue=0),
                 Coord.Cartesian(xmin=lb, ymin=lb, xmax=ub, ymax=ub),
                 Theme(panel_stroke="black"))
    return p
end

function get_data_metrics(data_df, metrics_df)
    concentrationBounds = percentile(data_df.Concentration, [2.5, 97.5])
    dataMetrics_df = combine(groupby(data_df, :exp_id),
                            :Viability => std => :viability_std,
                            :Concentration => minimum => :concentration_min,
                            :Concentration => maximum => :concentration_max)
    metrics_df = innerjoin(metrics_df, dataMetrics_df, on=:exp_id)
    return concentrationBounds, metrics_df
end

function ic50_std_plot(df, lb, ub, concentrationBounds)
    Gadfly.set_default_plot_size(3inch, 2inch)
    p = Gadfly.plot(df, x=:ic50, y=:viability_std, yintercept=20,
                 Geom.hline(), Geom.hexbin(xbincount=80, ybincount=80),
                 Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Coord.cartesian(xmin=lb, xmax=ub, ymin=0),
                 Theme(panel_stroke= "black"),
                 Guide.title("N="*string(length(unique(df.exp_id)))))

    push!(p, layer(xintercept=concentrationBounds, Geom.vline(color=["black"])))
    return p
end

function ic50_std_contour_plot(mtx, lb, ub, concentrationBounds)
    Gadfly.set_default_plot_size(3inch, 2inch)
    p = Gadfly.plot(z=mtx.density, x=mtx.x, y=mtx.y,
                    Geom.contour(levels=8),
                    Scale.color_continuous(colormap=scaleColor, minvalue=0),
                    Coord.cartesian(xmin=lb, xmax=ub, ymin=0),
                    Theme(panel_stroke= "black"))

    push!(p, layer(xintercept=concentrationBounds, Geom.vline(color=["black"])))
    return p
end

function get_prob_ic50(posterior_df, metrics_df)
    posterior_df = innerjoin(posterior_df, metrics_df[:,[:exp_id, :concentration_min, :concentration_max]], on=:exp_id)
posterior_df[:, :withinDose] = posterior_df.concentration_min .< posterior_df.ic50 .< posterior_df.concentration_max

metrics_df[:, :withinDoseML] = metrics_df.concentration_min .< metrics_df.ic50 .< metrics_df.concentration_max
metrics_df[:, :outsideDoseML] = .!metrics_df.withinDoseML

withinDose_prob = combine(groupby(posterior_df, :exp_id), :withinDose => sum => :withinCount,
                                                          :viability_std => unique => :viability_std)
withinDose_prob[:, :withinProb] = withinDose_prob.withinCount ./ 4000
withinDose_prob[:, :outsideProb] = 1. .- withinDose_prob.withinProb
withinDose_prob = innerjoin(withinDose_prob, metrics_df[:, [:exp_id, :convergence, :outsideDoseML]], on=:exp_id)
    return withinDose_prob, posterior_df
end

function boxplots_dose_plot(withinDose_prob)
    Gadfly.set_default_plot_size(2inch, 2inch)
    tmp = get_converged(withinDose_prob)
    pA = Gadfly.plot(withinDose_prob, x=:outsideDoseML, y=:outsideProb, 
                   Geom.boxplot(),)
    pA_conv = Gadfly.plot(tmp, x=:outsideDoseML, y=:outsideProb, 
                   Geom.boxplot(),)

    pB = Gadfly.plot(withinDose_prob, x=:convergence, y=:outsideProb, 
                    Geom.boxplot(),)

    pC = Gadfly.plot(withinDose_prob, x=:outsideDoseML, y=:viability_std, yintercept=[20],
                    Geom.hline(), Geom.boxplot())

    pC_conv = Gadfly.plot(tmp, x=:outsideDoseML, y=:viability_std, yintercept=[20],
                    Geom.hline(), Geom.boxplot())

    println("------> Count outsideDoseML : ", counter(withinDose_prob.outsideDoseML))
    println("------> Count outsideDoseML conv.: ", counter(tmp.outsideDoseML))
    println("------> Count convergence : ", counter(withinDose_prob.convergence))

    return pA, pA_conv, pB, pC, pC_conv
end

function get_posterior_diff(posterior_df, metrics_df)
    posteriorDiff_df = combine(groupby(posterior_df, :exp_id), 
                                    eff_metrics[1] => (x -> percentile(x, 97.5) - percentile(x, 2.5)) => :HDR_postDiff,
                                    eff_metrics[2] => (x -> percentile(x, 97.5) - percentile(x, 2.5)) => :LDR_postDiff,
                                    eff_metrics[3] => (x -> percentile(x, 97.5) - percentile(x, 2.5)) => :ic50_postDiff,
                                    eff_metrics[4] => (x -> percentile(x, 97.5) - percentile(x, 2.5)) => :slope_postDiff,
                                    eff_metrics[5] => (x -> percentile(x, 97.5) - percentile(x, 2.5)) => :aac_postDiff)

    metrics_df = innerjoin(metrics_df, posteriorDiff_df, on=:exp_id)
    return metrics_df
end

function ml_post_diff_plot(df, mt, lb, ub)
    Gadfly.set_default_plot_size(3inch, 2inch)
    p = Gadfly.plot(df, x=mt, y=Symbol(string(mt)*"_postDiff"), 
                    Geom.hexbin(xbincount=120, ybincount=100),
                    Scale.color_continuous(colormap=scaleColor, minvalue=1),
                    Coord.cartesian(xmin=lb, xmax=ub, ymin=0),
                    Theme(panel_stroke="black"))
    return p
end

function ml_post_diff_contour_plot(mtx, lb, ub)
    Gadfly.set_default_plot_size(3inch, 2inch)
    p = Gadfly.plot(z=mtx.density, x=mtx.x, y=mtx.y,
                 Geom.contour(levels=8),
                 Scale.color_continuous(colormap=scaleColor, minvalue=0),
                 Coord.Cartesian(xmin=lb, ymin=lb, xmax=ub, ymax=ub),
                 Theme(panel_stroke="black"))
    return p
end

function get_converged(df)
    tmp = filter(:convergence => x -> x == true, df)
    println("------> ", length(unique(tmp.exp_id)), " experiments converged out of ", length(unique(df.exp_id)))
    return tmp
end





#Gadfly.set_default_plot_size(3inch, 2inch)
#p6ya = Gadfly.plot(withinDose_prob, x=:viability_std, color=:convergence, 
#                   Geom.histogram(density=true),
#                   Coord.cartesian(xmin=0, xmax=60))
#draw(PDF(figure_prefix*"hist_conv_viabStd_ic50.pdf", 3inch, 2inch), p6ya)

#p6yb = Gadfly.plot(withinDose_prob, x=:viability_std, color=:outsideDoseML,
#                   Geom.histogram(density=true),
#                   Coord.cartesian(xmin=0, xmax=60))
#draw(PDF(figure_prefix*"hist_ml_viabStd_ic50.pdf", 3inch, 2inch), p6yb)


###### MAD vs. ML estimations ##########
#gCSImad_df = combine(groupby(gCSIposterior_df, :exp_id), 
#                        :HDR => mad => :HDR_mad,
#                        :LDR => mad => :LDR_mad,
#                        :ic50 => mad => :ic50_mad)
#filter(:exp_id => x -> x ∈ expIdSubset_list, gCSImad_df)

#gCSImetrics_df = innerjoin(gCSImetrics_df, gCSImad_df, on=:exp_id)


#### Example #4
## exp_id = HCC78_Lapatinib_11a, KP-3_Erastin_4b, KMM-1_Dabrafenib_11d
#e = "HCC78_Lapatinib_11a"
#Gadfly.set_default_plot_size(4inch, 2inch)
#viab_subset = filter(:exp_id => x -> x == e, data_df)
#p10 = Gadfly.plot(layer(viab_subset, x=:Concentration, y=:Viability, color=:exp_id, Geom.point()),
#                  Coord.cartesian(ymin=-10, ymax=110),
#                  Theme(panel_stroke="black"))

#subset_posterior = getBIDRAposterior("gCSI", [e])
#postCurves = getPosteriorCurves(subset_posterior, -4, 1)
#xDose = parse.(Float64, names(postCurves))

#med_curve = median.(eachcol(postCurves))
#upper_curve = percentile.(eachcol(postCurves), 97.5)
#lower_curve = percentile.(eachcol(postCurves), 2.5)
#push!(p10, layer(x=xDose, y=med_curve, color=[e], linestyle=[:dot], Geom.line()))
        
#mleFit = llogistic(filter(:exp_id => x -> x == e, gCSIml_df)[1, [:LDR, :HDR, :ic50, :slope]])
#yViability = mleFit.(xDose)
#push!(p10, layer(x=xDose, y=yViability, color=[e], Geom.line()))

#push!(p10, layer(x=xDose, ymin=lower_curve, ymax=upper_curve, color=[e], alpha=[0.3], Geom.ribbon()))
#draw(PDF(figure_prefix*"curveExample.pdf", 4inch, 2inch), p10)

#Gadfly.set_default_plot_size(3inch, 2inch)
#x = percentile(subset_posterior.HDR, [2.5, 97.5])
#ymin = [0.,0.]
#ymax = ymin .+ 0.04
#hdr = filter(:exp_id => x -> x == e, gCSIml_df)[:, :HDR]
#p11 = Gadfly.plot(subset_posterior, x=:HDR, Geom.histogram(density=true),
#                  layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon()),
#                  layer(xintercept=hdr, Geom.vline()),
#                  Theme(panel_stroke="black"))
#draw(PDF(figure_prefix*"posteriorExample_HDR.pdf", 3inch, 2inch), p11)

#Gadfly.set_default_plot_size(3inch, 2inch)
#x = percentile(subset_posterior.ic50, [2.5, 97.5])
#ymin = [0.,0.]
#ymax = ymin .+ 4
#ic50 = filter(:exp_id => x -> x == e, gCSIml_df)[:, :ic50]
#p12 = Gadfly.plot(subset_posterior, x=:ic50, Geom.histogram(density=true),
#                  layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon()),
#                  layer(xintercept=ic50, Geom.vline()),
#                  Theme(panel_stroke="black"))
#draw(PDF(figure_prefix*"posteriorExample_IC50.pdf", 3inch, 2inch), p12)