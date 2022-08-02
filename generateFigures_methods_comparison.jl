using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig
using Utils

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
    print("--> List of expId ")
    @time expId_list = getExpId_h5(dt);

    print("--> StrIndex ")
    @time si = Utils.StrIndex(expId_list);

    println("--> Raw data ")
    @time data_df = getRawData_h5(dt, false, si);

    println("--> ML estimates")
    @time ml_df = getMLestimates(dt, si);
    tmp = get_converged(ml_df);
    
    ###### Median vs. ML estimations ##########
    println("--> Posterior distributions")
    @time posterior_df = getPosterior_h5(dt, false, si);

    println("--> Add median and metrics")
    metrics_df = get_median(posterior_df, eff_metrics, ml_df);

    for em in eff_metrics
        println("---> ", em)
        println("------> Plotting hexbin")
        em_subset_df = metricZoom_subset(metrics_df, em, metrics_bounds[em][1], metrics_bounds[em][2])
        em_prior = getPriors[string(em)](N)
        p = medianML_hexbin_plot(em_subset_df, em, metrics_bounds[em][1], metrics_bounds[em][2], em_prior)

        if dt == "gCSI"
            exp_subset_df = filter(:exp_id => x -> si.id2str[x] ∈ expIdSubset_list, metrics_df)
            push!(p, layer(exp_subset_df, x=Symbol(String(em)*"_median"), y=em, Geom.point(), order=1))
        end

        draw(PDF(figure_prefix*dt*"_med_mle_"*string(em)*".pdf", 3inch, 2inch), p)

        println("------> Plotting contour")
        em_mtx = kde((em_subset_df[:, Symbol(String(em)*"_median")], em_subset_df[:, em]))
        p = medianML_contour_plot(em_mtx, metrics_bounds[em][1], metrics_bounds[em][2], em_prior)
        draw(PDF(figure_prefix*dt*"_med_mle_"*string(em)*"_contour.pdf", 3inch, 2inch), p)
    end
end