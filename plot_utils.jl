using DataFrames
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
    println("------> ", nrow(tmp), " experiments out of ", nrow(metrics_df), " for ", mt)
    return tmp
end

### Select posterior points that lay within a given range
function posteriorZoom_subset(posterior_df, mt, lb, ub)
    tmp = filter(mt => x -> lb ≤ x ≤ ub, posterior_df)
    println("------> ", length(unique(tmp.exp_id)), " experiments out of ", length(unique(posterior_df.exp_id)), " for ", mt)
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

function get_posterior_diff(posterior_dr, metrics_df)
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
