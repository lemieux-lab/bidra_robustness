using DataFrames, HDF5
using CSV
using Statistics, StatsBase
using Distributions
using Gadfly, StatsPlots
using Cairo, Fontconfig
using CairoMakie

include("utils.jl")

R = 100  # Number of repetitions for each condition
data_prefix = "_generated_data_reviews";
synthetic_fn = data_prefix*"/synthetic_data_gaussian_noise_results_$(string(R))_repeats.h5";

figure_prefix = "_generated_figures_reviews/synthetic_missing_doses/";
if !isdir(figure_prefix)
    mkpath(figure_prefix)
end

### Switches
plot_synthetic_data = true
plot_lsqfit_missing_data_effect = true
plot_bidra_missing_data_effect = true


### Initial true curve
n = 12
dose_range = Array(range(-3, stop=3, length=n));
plotting_dose_range = Array(range(-3, stop=3, length=1000));

efficiency_metrics = [100.0, 0.0, 0., 1.0];
true_curve_func = llogistic(efficiency_metrics);
true_response = true_curve_func.(dose_range);
true_response_plotting = true_curve_func.(plotting_dose_range);

### Synthetic data visualisation
if plot_synthetic_data
    dataset = h5open(synthetic_fn, "r")
    noise_groups = keys(dataset)

    for noise in noise_groups
        data_group = dataset[noise]["data"]
        subsets = keys(data_group)

        subset = "subset_12"
        subset_group = data_group[subset]
        repetitions = keys(subset_group)

        dose_response_data = reduce(vcat, [read(subset_group[rep]) for rep in repetitions])

        fig = CairoMakie.Figure(backgroundcolor="transparent", size=(400, 300));
        ax = CairoMakie.Axis(fig[1, 1], xlabel="Dose (log10)", ylabel="Viability (%)", title="Synthetic data with Gaussian noise (σ=$(split(noise, "_")[end]))")

        CairoMakie.violin!(ax, dose_response_data[:, 1], dose_response_data[:, 3], label="Noisy data", color=:black, width=0.4) ## Noisy dataset
        CairoMakie.scatter!(ax, dose_range, true_response, label="Original data", color=:red, markersize=8) ## Original dataset
        CairoMakie.lines!(ax, plotting_dose_range, true_response_plotting, label="True curve", color=:red, linewidth=1.5, linestyle=:dash, alpha=0.7) ## True curve

        fn = figure_prefix*"/complete_DR_$(string(noise))_$(string(R))_repetitions.pdf"
        CairoMakie.save(fn, fig)
    end
    close(dataset)
end


### Effect of missing data on lsqfit curves
SUBSET_ORDER = [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]
dataset = h5open(synthetic_fn, "r")
noise_groups = keys(dataset)

for noise in noise_groups
    estimate_group = dataset[noise]["estimates"]
    subsets = keys(estimate_group)

    if plot_lsqfit_missing_data_effect
        fig = CairoMakie.Figure(backgroundcolor="transparent", size=(450*11, 350*1));

        for subset in subsets
            subset_int = parse(Int, split(subset, "_")[end])
            idx = findfirst(val -> val == subset_int, SUBSET_ORDER)
            
            # Skip if subset not in order list
            if isnothing(idx)
                continue
            end
            
            col = idx * 2
            subset_group = estimate_group[subset]
            repetitions = keys(subset_group)

            ### LsqFit results
            row = 2
            ax_upper = CairoMakie.Axis(fig[row-1, col], title="Synthetic data with Gaussian noise (σ=$(split(noise, "_")[end]))\nSubset with $(subset_int) data points",
                                    xticksvisible=false, yticksvisible=false, xticklabelsvisible=false, yticklabelsvisible=false)
            ic50 = [read(estimate_group[subset][rep])[3] for rep in repetitions ]
            ic50_clamp = clamp.(ic50, -3.5, 3.5)
            CairoMakie.hist!(ax_upper, ic50_clamp, color=:black, bins=25)
            CairoMakie.vlines!(ax_upper, [0.], color=:red, linestyle=:dash) ## True ic50

            ax_right = CairoMakie.Axis(fig[row, col+1], 
                                    xticksvisible=false, yticksvisible=false, xticklabelsvisible=false, yticklabelsvisible=false)
            hdr = [read(estimate_group[subset][rep])[1] for rep in repetitions ]
            ldr = [read(estimate_group[subset][rep])[2] for rep in repetitions ]
            hdr_clamp = clamp.(hdr, -50, 150)
            ldr_clamp = clamp.(ldr, -50, 150)
            CairoMakie.hist!(ax_right, hdr_clamp, color=:black, direction=:x, bins=25)
            CairoMakie.hist!(ax_right, ldr_clamp, color=:gray, direction=:x, bins=25)
            CairoMakie.hlines!(ax_right, [0., 100.], color=:red, linestyle=:dash) ## True ic50

            ax_main = CairoMakie.Axis(fig[row, col], xlabel="Dose (log10)", ylabel="Viability (%)")
            for rep in repetitions
                rep_estimate = read(estimate_group[subset][rep])
                rep_curve_func = llogistic(rep_estimate)
                rep_response = rep_curve_func.(plotting_dose_range)
                rep_response = clamp.(rep_response, -50, 150)  # Clamp to visible range

                CairoMakie.lines!(ax_main, plotting_dose_range, rep_response, color=:gray, linewidth=1.5, alpha=0.4) ## LsqFit curve
            end

            ### Add data points
            data = dataset[noise]["data"][subset]
            dose_response_data = reduce(vcat, [read(data[rep]) for rep in repetitions])
            CairoMakie.violin!(ax_main, dose_response_data[:, 1], dose_response_data[:, 3], label="Noisy data", color=:black, width=0.4) ## Noisy dataset

            if subset_int < 12 && subset_int >= 1 && subset_int+1 <= length(dose_range)
                CairoMakie.scatter!(ax_main, dose_range[subset_int+1:end], true_response[subset_int+1:end], label="Missing data", color=:red, markersize=8) ## Missing data points
            end   

            CairoMakie.limits!(ax_main, -3.5, 3.5, -50, 150)
            CairoMakie.xlims!(ax_upper, -3.5, 3.5)
            CairoMakie.ylims!(ax_right, -50, 150)

            CairoMakie.rowsize!(fig.layout, row-1, 50)
            CairoMakie.rowgap!(fig.layout, row-1, 0)

            CairoMakie.colsize!(fig.layout, col+1, 50)
            CairoMakie.colgap!(fig.layout, col, 0)
        end

        fn = figure_prefix*"/effect_of_missing_data_on_LM_$(string(noise))_$(string(R))_repetitions.pdf"
        CairoMakie.save(fn, fig)

    end


    ## BiDRA results
    function get_posterior_distribution(group, col)
        repetitions = keys(group)
        posterior_values = [read(group[rep])[:, col] for rep in repetitions ]
        return vcat(posterior_values...)
    end

    if plot_bidra_missing_data_effect
        for posterior_name in ["posterior", "no_prior_posterior"]
            posterior_group = dataset[noise][posterior_name]
            subsets = keys(estimate_group)

            fig = CairoMakie.Figure(backgroundcolor="transparent", size=(450*11, 350*1));

            for subset in subsets
                subset_int = parse(Int, split(subset, "_")[end])
                idx = findfirst(val -> val == subset_int, SUBSET_ORDER)
                
                # Skip if subset not in order list
                if isnothing(idx)
                    continue
                end
                
                col = idx * 2
                subset_group = posterior_group[subset]
                repetitions = keys(subset_group)

                row = 2
                ax_upper = CairoMakie.Axis(fig[row-1, col], title="Synthetic data with Gaussian noise (σ=$(split(noise, "_")[end]))\nSubset with $(subset_int) data points",
                                        xticksvisible=false, yticksvisible=false, xticklabelsvisible=false, yticklabelsvisible=false)
                
                ic50 = get_posterior_distribution(subset_group, 3)
                ic50_clamp = clamp.(ic50, -3.5, 3.5)
                CairoMakie.hist!(ax_upper, ic50_clamp, color=:black, bins=25)
                CairoMakie.vlines!(ax_upper, [0.], color=:red, linestyle=:dash) ## True ic50

                ax_right = CairoMakie.Axis(fig[row, col+1], 
                                            xticksvisible=false, yticksvisible=false, xticklabelsvisible=false, yticklabelsvisible=false)
                hdr = get_posterior_distribution(subset_group, 1)
                ldr = get_posterior_distribution(subset_group, 2)
                hdr_clamp = clamp.(hdr, -50, 150)
                ldr_clamp = clamp.(ldr, -50, 150)
                CairoMakie.hist!(ax_right, hdr_clamp, color=:black, direction=:x, bins=25)
                CairoMakie.hist!(ax_right, ldr_clamp, color=:gray, direction=:x, bins=25)
                CairoMakie.hlines!(ax_right, [0., 100.], color=:red, linestyle=:dash) ## True extremes

                ax_main = CairoMakie.Axis(fig[row, col], xlabel="Dose (log10)", ylabel="Viability (%)")
                for rep in repetitions
                    rep_posterior = rand(eachrow(read(subset_group[rep])))
                    rep_curve_func = llogistic(rep_posterior[1:4])
                    rep_response = rep_curve_func.(plotting_dose_range)
                    rep_response = clamp.(rep_response, -50, 150)  # Clamp to visible range

                    CairoMakie.lines!(ax_main, plotting_dose_range, rep_response, color=:gray, linewidth=1.5, alpha=0.4) ## LsqFit curve
                end

                ### Add data points
                data = dataset[noise]["data"][subset]
                dose_response_data = reduce(vcat, [read(data[rep]) for rep in repetitions])
                CairoMakie.violin!(ax_main, dose_response_data[:, 1], dose_response_data[:, 3], label="Noisy data", color=:black, width=0.4) ## Noisy dataset

                if subset_int < 12 && subset_int >= 1 && subset_int+1 <= length(dose_range)
                    CairoMakie.scatter!(ax_main, dose_range[subset_int+1:end], true_response[subset_int+1:end], label="Missing data", color=:red, markersize=8) ## Missing data points
                end   

                CairoMakie.limits!(ax_main, -3.5, 3.5, -50, 150)
                CairoMakie.xlims!(ax_upper, -3.5, 3.5)
                CairoMakie.ylims!(ax_right, -50, 150)

                CairoMakie.rowsize!(fig.layout, row-1, 50)
                CairoMakie.rowgap!(fig.layout, row-1, 0)

                CairoMakie.colsize!(fig.layout, col+1, 50)
                CairoMakie.colgap!(fig.layout, col, 0)
            end

            fn = figure_prefix*"/effect_of_missing_data_on_BiDRA_$(posterior_name)_$(string(noise))_$(string(R))_repetitions.pdf"
            CairoMakie.save(fn, fig)
        end
    end
end

close(dataset)
