using DataFrames, HDF5
using CSV
using Statistics, StatsBase
using Distributions
using Gadfly, StatsPlots
using Cairo, Fontconfig
using CairoMakie

include("utils.jl")
#include("reviews_generateFigure_synthetic_data.jl")

function get_posterior_distribution(group, col)
    repetitions = keys(group)
    posterior_values = [read(group[rep])[:, col] for rep in repetitions]
    return vcat(posterior_values...)
end

R = 100  # Number of repetitions for each condition
data_prefix = "_generated_data_reviews";
synthetic_fn = data_prefix*"/synthetic_data_gaussian_noise_results_$(string(R))_repeats.h5";

figure_prefix = "_generated_figures_reviews/compare_metrics_to_true/";
if !isdir(figure_prefix)
    mkpath(figure_prefix)
end

### Initial true curve
n = 12
dose_range = Array(range(-3, stop=3, length=n));
plotting_dose_range = Array(range(-3, stop=3, length=1000));

efficiency_metrics = Dict("IC50" => Dict("true" => 0.0, "lower_bound" => -3.5, "upper_bound" => 3.5, "idx" => 3),
                          "LDR" => Dict("true" => 100.0, "lower_bound" => -50.0, "upper_bound" => 150.0, "idx" => 1),
                          "Slope" => Dict("true" => 1.0, "lower_bound" => 0.0, "upper_bound" => 5.0, "idx" => 4),
                          "HDR" => Dict("true" => 0.0, "lower_bound" => -50.0, "upper_bound" => 150.0, "idx" => 2));
true_efficiency_metrics = [100.0, 0.0, 0., 1.0];
true_curve_func = llogistic(true_efficiency_metrics);
true_response = true_curve_func.(dose_range);
true_response_plotting = true_curve_func.(plotting_dose_range);

### Compare efficiency metrics to true values
SUBSET_ORDER = [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]
dataset = h5open(synthetic_fn, "r")
noise_groups = keys(dataset)

for noise in noise_groups
    estimate_group = dataset[noise]["estimates"]
    bidra_group = dataset[noise]["posterior"]
    no_prior_group = dataset[noise]["no_prior_posterior"]
    subsets = keys(estimate_group)

    for metric_val in keys(efficiency_metrics)
        metric_info = efficiency_metrics[metric_val]
        println("Noise: $(noise), Metric: $(metric_val), True value: $(metric_info["true"])")

        

        for subset in subsets
            subset_int = parse(Int, split(subset, "_")[end])
            idx = findfirst(val -> val == subset_int, SUBSET_ORDER)

            # Skip if subset not in order list
            if isnothing(idx)
                continue
            end
            
            col = idx
            #metric_info = efficiency_metrics["IC50"]

            fig = CairoMakie.Figure(backgroundcolor="transparent", size=(450, 350));
            ax = CairoMakie.Axis(fig[1, 1], xlabel="Delta $(metric_val)", ylabel="Method", title="Subset $(subset_int) doses")

            ### Delta between estimated and true metric values
            subset_group = estimate_group[subset]
            repetitions = keys(subset_group)

            metrics = [read(estimate_group[subset][rep])[metric_info["idx"]] for rep in repetitions ]
            metric_clamp = clamp.(metrics, metric_info["lower_bound"], metric_info["upper_bound"])
            metric_delta = metric_clamp .- metric_info["true"]

            CairoMakie.rainclouds!(ax, repeat([3], length(metric_delta)), metric_delta, 
                                   clouds=hist, markersize=0.0, orientation=:horizontal)

            ### Delta between inferred and true metric values
            subset_group = bidra_group[subset]

            metrics = get_posterior_distribution(subset_group, metric_info["idx"])
            metric_clamp = clamp.(metrics, metric_info["lower_bound"], metric_info["upper_bound"])
            metric_delta = metric_clamp .- metric_info["true"]

            CairoMakie.rainclouds!(ax, repeat([2], length(metric_delta)), metric_delta, 
                                   clouds=hist, markersize=0.0, orientation=:horizontal)

            ### Delta between inferred and true metric values (no prior)
            subset_group = no_prior_group[subset]

            metrics = get_posterior_distribution(subset_group, metric_info["idx"])
            metric_clamp = clamp.(metrics, metric_info["lower_bound"], metric_info["upper_bound"])
            metric_delta = metric_clamp .- metric_info["true"]

            CairoMakie.rainclouds!(ax, repeat([1], length(metric_delta)), metric_delta, 
                                   clouds=hist, markersize=0.0, orientation=:horizontal)

            #fn = figure_prefix*"/compare_true_and_predicted_$(string(metric_val))_with_$(string(subset_int))_response_$(string(noise))_$(string(R))_repetitions.png"
            #CairoMakie.save(fn, fig)

            fn = figure_prefix*"/compare_true_and_predicted_$(string(metric_val))_with_$(string(subset_int))_response_$(string(noise))_$(string(R))_repetitions.pdf"
            CairoMakie.save(fn, fig)
        end
    end

end
