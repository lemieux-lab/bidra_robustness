using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")
include("compound_characterization/MCMCmodels.jl")

data_prefix = "_generated_data_reviews/synthetic_data";

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

function do_curve_fit_with_LsqFit(data_df::DataFrame)
    @. model(x, p) = p[2] + ((p[1] - p[2]) / (1 + 10^(p[4] * (x - p[3]))))
    p₀ = [100.,0.,0.,1]

    fit = curve_fit(model, data_df.x, data_df.y_noisy, p₀)

    bestFitParam = fit.param
    rmse = sqrt(sum(fit.resid .^ 2) / length(fit.resid))
    convergence = fit.converged
    return bestFitParam, rmse, convergence
end

for σ in [0.1, 1., 5., 10.]
    sigma_path = data_prefix*"/synthetic_data_gaussian_noise_$(string(σ))"
    if !isdir(sigma_path)
        mkdir(sigma_path)
    end

    n = 12
    dose_range = Array(range(-3, stop=3, length=n));
    efficiency_metrics = [100.0, 0.0, 0., 1.0];
    synthetic_df = make_dataset_with_gaussian_noise(dose_range, efficiency_metrics, σ);

    fn = sigma_path*"/synthetic_data_gaussian_noise_$(string(σ)).csv";
    CSV.write(fn, synthetic_df);
    
    posterior_df = do_inference_with_BiDRA(synthetic_df);
    CSV.write(sigma_path*"/synthetic_data_gaussian_noise_$(string(σ))_posterior.csv", posterior_df);

    estimates, rmse, convergence_status = do_curve_fit_with_LsqFit(synthetic_df);
    estimate_df = DataFrame(LDR=estimates[1], HDR=estimates[2], ic50=estimates[3],
                        slope=estimates[4], rmse=rmse, convergence=convergence_status);
    CSV.write(sigma_path*"/synthetic_data_gaussian_noise_$(string(σ))_estimates.csv", estimate_df);

    for r in 1:100 
        println("σ=$(string(σ)), repetition $r/100")
    end

    for i in n-1:-1:2
        subset_path = sigma_path*"/subset_$(string(i))_response"
        if !isdir(subset_path)
            mkdir(subset_path)
        end

        synthetic_df_subset = synthetic_df[1:i, :]
        fn = subset_path*"/synthetic_data_gaussian_noise_$(string(σ))_subset_$(string(i))_datapoints.csv";
        CSV.write(fn, synthetic_df_subset);

        posterior_fn = subset_path*"/synthetic_data_gaussian_noise_$(string(σ))_datapoints_posterior.h5";
        lsqfit_fn = subset_path*"/synthetic_data_gaussian_noise_$(string(σ))_datapoints_estimates.h5";

        for r in 1:100 
            println("σ=$(string(σ)), $(string(i)) datapoints, repetition $r/100")
            posterior_df = do_inference_with_BiDRA(synthetic_df_subset);

            estimates, rmse, convergence_status = do_curve_fit_with_LsqFit(synthetic_df_subset);
            estimate_df = DataFrame(LDR=estimates[1], HDR=estimates[2], ic50=estimates[3],
                            slope=estimates[4], rmse=rmse, convergence=convergence_status);

            h5open(posterior_fn, "cw") do file
                g = create_group(file, "repetition_$(r)")
                g["posterior"] = convert(Array{Float64}, Array(posterior_df))
                g["median"] = median.(eachcol(posterior_df[:,[:HDR, :LDR, :ic50, :slope, :σ]]))
            end

            h5open(lsqfit_fn, "cw") do file
                g = create_group(file, "repetition_$(r)")
                g["estimates"] = convert(Array{Float64}, Array(estimate_df))
            end

        end

    end
end