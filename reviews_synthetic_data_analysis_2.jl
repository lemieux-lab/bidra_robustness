using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")
include("compound_characterization/MCMCmodels.jl")

data_prefix = "_generated_data_reviews";

function make_dataset_with_gaussian_noise(synthetic_dose::Array, efficiency_metrics::Array, noise_std::Float64)
    ### Parameters: LDR, HDR, ic50, slope
    curve_func = llogistic(efficiency_metrics)
    synthetic_response = curve_func.(synthetic_dose)
    synthetic_response_with_noise = synthetic_response .+ rand(Normal(0, noise_std), length(synthetic_response));

    return DataFrame(x=synthetic_dose, y=synthetic_response, y_noisy=synthetic_response_with_noise);
end

function do_inference_with_BiDRA(data_df::DataFrame)
    nChain = 4
    nIte = 1000
    nAdapt = 1000
    δ = 0.65

    bidra_model = BIDRA(data_df.x, data_df.y_noisy)
    bidra_sampler = NUTS(nAdapt, δ)
    bidra_chains = sample(bidra_model, bidra_sampler, MCMCThreads(), nIte, nChain)

    posterior_df = DataFrame(bidra_chains)[:,[:LDR, :HDR, :ic50, :slope, :σ, :chain]]
    return posterior_df
end

function do_inference_with_BiDRA_no_prior(data_df::DataFrame)
    nChain = 4
    nIte = 1000
    nAdapt = 1000
    δ = 0.65

    p₀ = (HDR=0.0, LDR=100.0, ic50=0.0, slope=1.0, σ=1.0)

    bidra_model = BIDRA_no_prior(data_df.x, data_df.y_noisy)
    bidra_sampler = NUTS(nAdapt, δ)
    bidra_chains = sample(bidra_model, bidra_sampler, MCMCThreads(), nIte, nChain, initial_params=fill(Turing.InitFromParams(p₀), nChain))

    posterior_df = DataFrame(bidra_chains)[:,[:LDR, :HDR, :ic50, :slope, :σ, :chain]]
    return posterior_df
end

function do_curve_fit_with_LsqFit(data_df::DataFrame)
    @. model(x, p) = p[2] + ((p[1] - p[2]) / (1 + 10^(p[4] * (x - p[3]))))
    p₀ = [100.,0.,0.,1]

    fit = curve_fit(model, data_df.x, data_df.y_noisy, p₀)

    estimates = fit.param
    rmse = sqrt(sum(fit.resid .^ 2) / length(fit.resid))
    convergence = fit.converged

    estimate_df = DataFrame(LDR=estimates[1], HDR=estimates[2], ic50=estimates[3],
                            slope=estimates[4], rmse=rmse, convergence=convergence);

    return estimate_df
end

n = 12
dose_range = Array(range(-3, stop=3, length=n));
efficiency_metrics = [100.0, 0.0, 0., 1.0];

R = 100  # Number of repetitions for each condition

synthetic_fn = data_prefix*"/synthetic_data_gaussian_noise_results_$(string(R))_repeats.h5";
if isfile(synthetic_fn)
    rm(synthetic_fn)
end

h5open(synthetic_fn, "cw") do file
    for σ in [0.1, 1., 5., 10.]
        sigma_group = create_group(file, "sigma_$(string(σ))")

        data_group = create_group(sigma_group, "data")
        posterior_group = create_group(sigma_group, "posterior")
        no_prior_posterior_group = create_group(sigma_group, "no_prior_posterior")
        estimate_group = create_group(sigma_group, "estimates")

        for r in 1:R
            println("σ=$(string(σ)), repetition $r/$R")
            synthetic_df = make_dataset_with_gaussian_noise(dose_range, efficiency_metrics, σ);

            for i in n:-1:2
                subset_df = synthetic_df[1:i, :]
                posterior_df = do_inference_with_BiDRA(subset_df);
                no_prior_posterior_df = do_inference_with_BiDRA_no_prior(subset_df);
                estimate_df = do_curve_fit_with_LsqFit(subset_df);

                if r == 1
                    data_subset = create_group(data_group, "subset_$(string(i))")
                    posterior_subset = create_group(posterior_group, "subset_$(string(i))")
                    no_prior_posterior_subset = create_group(no_prior_posterior_group, "subset_$(string(i))")
                    estimate_subset = create_group(estimate_group, "subset_$(string(i))")
                end

                data_group["subset_$(string(i))"]["repetition_$(string(r))"] = convert(Array{Float64}, Array(subset_df))
                posterior_group["subset_$(string(i))"]["repetition_$(string(r))"] = convert(Array{Float64}, Array(posterior_df))
                no_prior_posterior_group["subset_$(string(i))"]["repetition_$(string(r))"] = convert(Array{Float64}, Array(no_prior_posterior_df))
                estimate_group["subset_$(string(i))"]["repetition_$(string(r))"] = convert(Array{Float64}, Array(estimate_df))
            end
        end
    end
end
