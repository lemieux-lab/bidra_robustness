push!(LOAD_PATH, "Utils/")
using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig
using Utils

include("utils.jl")
include("plot_utils.jl")

###### Global Var ####
figure_prefix = "_generated_figures/robustness/";

datasets = ["gray", "gCSI", "ctrpv2"];
eff_metrics = [:HDR, :LDR, :ic50, :slope, :aac];
metrics_bounds = Dict(:HDR=>[-50,150], :LDR=>[70,150], :ic50=>[-10,10], :slope=>[0,10], :aac=>[0,100]);
expIdSubset_list = ["NCI-H1648_AZ-628_8h", "Calu-1_PF-4708671_6b", "RERF-LC-MS_Gemcitabine_4b", "HCC78_Lapatinib_11a"];
scaleColor = Scale.lab_gradient("gray95","black");

for i in 1:length(datasets)
    dt = datasets[i];
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

    println("--> Posterior distributions")
    @time posterior_df = getPosterior_h5(dt, false, si);

    println("--> Add median and metrics")
    metrics_df = get_median(posterior_df, eff_metrics, ml_df);

    ###### IC50 vs. experimental concentration ##########
    println("2. Looking @ IC50 and experimental dose range")
    concentrationBounds, metrics_df = get_data_metrics(data_df, metrics_df);
    posterior_df = innerjoin(posterior_df, metrics_df[:,[:exp_id, :viability_std]], on=:exp_id);

    ### compare ic50 estimation to std
    println("---> Plotting IC50 est. vs. std")
    subset_df = metricZoom_subset(metrics_df, :ic50, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2]);

    println("------> Size of subset: ", size(subset_df))
    p = ic50_std_plot(subset_df, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2], concentrationBounds);
    if dt == "gCSI"
        exp_subset_df = filter(:exp_id => x -> si.id2str[x] ∈ expIdSubset_list, subset_df)
        push!(p, layer(exp_subset_df, x=:ic50, y=:viability_std, Geom.point(), order=1))
    end
    draw(PDF(figure_prefix*dt*"_std_mle_ic50.pdf", 3inch, 2inch), p);

    println("---> Plotting IC50 posterior vs. std")
    subset_df = metricZoom_subset(posterior_df, :ic50, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2]);
    p = ic50_std_plot(posterior_df, metrics_bounds[:ic50][1], metrics_bounds[:ic50][2], concentrationBounds);
    if dt == "gCSI"
        exp_subset_df = filter(:exp_id => x -> si.id2str[x] ∈ expIdSubset_list, subset_df)
        tmp_interval = combine(groupby(exp_subset_df, :exp_id), :ic50 => (x -> percentile(x, [2.5, 97.5])) => :interval)
        tmp_std = combine(groupby(exp_subset_df, :exp_id), :viability_std => unique => :viability_std)
        tmp = innerjoin(tmp_interval, tmp_std, on=:exp_id)
        push!(p, layer(tmp, x=:interval, y=:viability_std, group=:exp_id, Geom.line(), order=1))
    end
    draw(PDF(figure_prefix*dt*"_std_posterior_ic50.pdf", 3inch, 2inch), p);

    println("---> Calculating Prob. of IC50 being outside of experimental dose range")
    withinDose_prob, posterior_df = get_prob_ic50(posterior_df, metrics_df)

    println("---> Plotting prob. vs. std")
    Gadfly.set_default_plot_size(5inch, 2inch)
    p = Gadfly.plot(withinDose_prob, x=:outsideProb, y=:viability_std, yintercept=[20],
                    Geom.hline, Geom.hexbin(xbincount=120, ybincount=90),
                    Scale.color_continuous(colormap=scaleColor, minvalue=1),
                    Coord.cartesian(xmin=0, ymin=0, xmax=1),
                    Theme(panel_stroke="black"));

    if dt == "gCSI"
        exp_subset_df = filter(:exp_id => x -> si.id2str[x] ∈ expIdSubset_list, withinDose_prob)
        push!(p, layer(exp_subset_df, x=:outsideProb, y=:viability_std, Geom.point()))
    end
    draw(PDF(figure_prefix*dt*"_std_prob_ic50.pdf", 5inch, 2inch), p)

    println("---> Plotting series of boxplot")
    pA, pA_conv, pB, pC, pC_conv = boxplots_dose_plot(withinDose_prob);
    draw(PDF(figure_prefix*dt*"_boxplot_std_ml_ic50.pdf", 3inch, 3inch), pC)
    draw(PDF(figure_prefix*dt*"_boxplot_std_ml_ic50_conv.pdf", 3inch, 3inch), pC_conv)

    ###### Metrics vs. posterior ##########
    println("4. Looking @ HDR ML est. vs. posterior ΔHDR")
    metrics_df = get_posterior_diff(posterior_df, metrics_df)

    em = :HDR
    println("---> ", em)
    println("------> Plotting hexbin")
    em_subset_df = metricZoom_subset(metrics_df, em, metrics_bounds[em][1], metrics_bounds[em][2])
    p = ml_post_diff_plot(em_subset_df, em, metrics_bounds[em][1], metrics_bounds[em][2])

    if dt == "gCSI"
        exp_subset_df = filter(:exp_id => x -> si.id2str[x] ∈ expIdSubset_list, metrics_df)
        push!(p, layer(exp_subset_df, x=em, y=Symbol(String(em)*"_postDiff"), Geom.point(), order=1))
    end
    draw(PDF(figure_prefix*dt*"_postDiff_ml_"*string(em)*".pdf", 3inch, 2inch), p)

    ### compare HDR estimation to std
    Gadfly.set_default_plot_size(3inch, 2inch)
    p = Gadfly.plot(em_subset_df, x=:HDR, y=:viability_std, yintercept=[20],
                    Geom.hline, Geom.hexbin(xbincount=120, ybincount=90),
                    Scale.color_continuous(colormap=scaleColor, minvalue=1),
                    Coord.cartesian(xmin=metrics_bounds[em][1], xmax=metrics_bounds[em][2]),
                    Guide.title("n="*string(nrow(em_subset_df))),
                    Theme(panel_stroke="black"));

    if dt == "gCSI"
        exp_subset_df = filter(:exp_id => x -> si.id2str[x] ∈ expIdSubset_list, em_subset_df)
        push!(p, layer(exp_subset_df, x=:HDR, y=:viability_std, Geom.point()))
    end
    draw(PDF(figure_prefix*dt*"_std_ml_HDR.pdf", 3inch, 2inch), p)

    ### compare HDR posterior to std
    em_subset_df = metricZoom_subset(metrics_df, Symbol(string(em)*"_postDiff"), metrics_bounds[em][1], metrics_bounds[em][2])
    Gadfly.set_default_plot_size(3inch, 2inch)
    p = Gadfly.plot(em_subset_df, x=Symbol(string(em)*"_postDiff"), y=:viability_std, yintercept=[20],
                    Geom.hline, Geom.hexbin(xbincount=120, ybincount=90),
                    Scale.color_continuous(colormap=scaleColor, minvalue=1),
                    Coord.cartesian(xmin=metrics_bounds[em][1], xmax=metrics_bounds[em][2]),
                    Guide.title("n="*string(nrow(em_subset_df))),
                    Theme(panel_stroke="black"));

    if dt == "gCSI"
        exp_subset_df = filter(:exp_id => x -> si.id2str[x] ∈ expIdSubset_list, em_subset_df)
        push!(p, layer(exp_subset_df, x=Symbol(string(em)*"_postDiff"), y=:viability_std, Geom.point()))
    end
    draw(PDF(figure_prefix*dt*"_std_posterior_HDR.pdf", 3inch, 2inch), p)
end
