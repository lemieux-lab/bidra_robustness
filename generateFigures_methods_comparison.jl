using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")
include("plot_utils.jl")

###### Global Var ####
figure_prefix = "_generated_figures/methods_comparison/"
datasets = ["gray", "gCSI", "ctrpv2"]  
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
for i in 1:length(datasets)
    dt = datasets[i]
    println("**********************", dt, "**********************")
    println("1. Get data and metrics")
    @time data_df = getRawData_h5(dt)
    @time expId_list = getExpId_h5(dt)#unique(data_df.exp_id)

    @time ml_df = getMLestimates([dt], false, missing)
    tmp = get_converged(ml_df)

    @time posterior_df = getPosterior_h5(dt)#getBIDRAposterior(dt, expId_list)
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
    mtx = kde((subset_df.ic50, subset_df.viability_std))
    p = ic50_std_contour_plot(mtx, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2], concentrationBounds)
    draw(PDF(figure_prefix*dt*"_std_mle_ic50_contour.pdf", 3inch, 2inch), p)

    println("---> Plotting IC50 posterior vs. std")
    subset_df = metricZoom_subset(posterior_df, :ic50, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2])
    p = ic50_std_plot(subset_df, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2], concentrationBounds)
    draw(PDF(figure_prefix*dt*"_std_posterior_ic50.pdf", 3inch, 2inch), p)

    println("---> Plotting IC50 posterior vs. std contour")
    mtx = kde((subset_df.ic50, subset_df.viability_std))
    p = ic50_std_contour_plot(mtx, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2], concentrationBounds)
    draw(PDF(figure_prefix*dt*"_std_posterior_ic50_contour.pdf", 3inch, 2inch), p)

    println("---> Calculating Prob. of IC50 being outside of experimental dose range")
    withinDose_prob, posterior_df = get_prob_ic50(posterior_df, metrics_df)

    println("---> Plotting prob. vs. std")
    Gadfly.set_default_plot_size(5inch, 2inch)
    p = Gadfly.plot(withinDose_prob, x=:outsideProb, y=:viability_std, yintercept=[20],
                    Geom.hline, Geom.hexbin(xbincount=120, ybincount=90),
                    Scale.color_continuous(colormap=scaleColor, minvalue=1),
                    Coord.cartesian(xmin=0, ymin=0, xmax=1),
                    Theme(panel_stroke="black"))

    if dt == datasets[1]
        exp_subset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, withinDose_prob)
        push!(p, layer(exp_subset_df, x=:outsideProb, y=:viability_std, Geom.point()))
    end
    draw(PDF(figure_prefix*dt*"_std_prob_ic50.pdf", 5inch, 2inch), p)

    println("---> Plotting prob. vs. std contour")
    mtx = kde((withinDose_prob.outsideProb, withinDose_prob.viability_std))
    p = Gadfly.plot(z=mtx.density, x=mtx.x, y=mtx.y,
                    Geom.contour(levels=8),
                    Scale.color_continuous(colormap=scaleColor, minvalue=0),
                    Coord.cartesian(xmin=0, ymin=0, xmax=1),
                    Theme(panel_stroke="black"))
    draw(PDF(figure_prefix*dt*"_std_prob_ic50_contour.pdf", 5inch, 2inch), p)


    println("---> Plotting series of boxplot")
    pA, pA_conv, pB, pC, pC_conv = boxplots_dose_plot(withinDose_prob)
    draw(PDF(figure_prefix*dt*"_boxplot_ml_prob_ic50.pdf", 3inch, 3inch), pA)
    draw(PDF(figure_prefix*dt*"_boxplot_ml_prob_ic50_conv.pdf", 3inch, 3inch), pA_conv)
    draw(PDF(figure_prefix*dt*"boxplot_conv_prob_ic50.pdf", 3inch, 3inch), pB)
    draw(PDF(figure_prefix*dt*"boxplot_std_ml_ic50.pdf", 3inch, 3inch), pC)
    draw(PDF(figure_prefix*dt*"boxplot_std_ml_ic50_conv.pdf", 3inch, 3inch), pC_conv)


    ###### Metrics vs. posterior ##########
    println("4. Looking @ ML est. vs. posterior Δ")
    metrics_df = get_posterior_diff(posterior_df, metrics_df)

    for em in eff_metrics
        println("---> ", em)
        println("------> Plotting hexbin")
        em_subset_df = metricZoom_subset(metrics_df, em, metrics_bounds[em][1], metrics_bounds[em][2])
        p = ml_post_diff_plot(em_subset_df, em, metrics_bounds[em][1], metrics_bounds[em][2])

        if dt == datasets[1]
            exp_subset_df = filter(:exp_id => x -> x ∈ expIdSubset_list, metrics_df)
            push!(p, layer(exp_subset_df, x=Symbol(String(em)*"_postDiff"), y=em, Geom.point()))
        end

        draw(PDF(figure_prefix*dt*"_postDiff_ml_"*string(em)*".pdf", 3inch, 2inch), p)

        println("------> Plotting contour")
        em_mtx = kde((em_subset_df[:, Symbol(String(em)*"_postDiff")], em_subset_df[:, em]))
        p = ml_post_diff_contour_plot(em_mtx, metrics_bounds[em][1], metrics_bounds[em][2])
        draw(PDF(figure_prefix*dt*"_postDiff_ml_"*string(em)*"_contour.pdf", 3inch, 2inch), p)
    end
end