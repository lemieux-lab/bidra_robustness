using CSV, DataFrames
using Distributions, Statistics, StatsBase
using LsqFit
using Logging

include("compound_characterization/MCMCmodels.jl")

data_prefix = "_generated_data_reviews";

const N_CHAIN = 4
const N_ITE = 1000000
const N_ADAPT = 1000
const TARGET_ACCEPT = 0.65
const LSQFIT_P0 = [100.0, 0.0, 0.0, 1.0]

@. lsqfit_model(x, p) = p[2] + ((p[1] - p[2]) / (1 + 10^(p[4] * (x - p[3]))))

function make_dataset_with_gaussian_noise(synthetic_dose::Array, efficiency_metrics::Array, noise_std::Float64)
    ### Parameters: LDR, HDR, ic50, slope
    curve_func = llogistic(efficiency_metrics)
    synthetic_response = curve_func.(synthetic_dose)
    synthetic_response_with_noise = synthetic_response .+ rand(Normal(0, noise_std), length(synthetic_response));

    return DataFrame(x=synthetic_dose, y=synthetic_response, y_noisy=synthetic_response_with_noise);
end

function do_inference_with_BiDRA(xs, ys)
    bidra_model = BIDRA(xs, ys)
    bidra_sampler = NUTS(N_ADAPT, TARGET_ACCEPT)

    # Silence sampler informational logs (e.g. initial step size) during benchmarks.
    return @elapsed with_logger(NullLogger()) do
        sample(bidra_model, bidra_sampler, MCMCThreads(), N_ITE, N_CHAIN; progress=false, verbose=false)
    end

    #posterior_df = DataFrame(bidra_chains)[:,[:LDR, :HDR, :ic50, :slope, :σ, :chain]]
    #return posterior_df
end

function do_curve_fit_with_LsqFit(xs, ys)
    return @elapsed curve_fit(lsqfit_model, xs, ys, LSQFIT_P0)

    #estimates = fit.param
    #rmse = sqrt(sum(fit.resid .^ 2) / length(fit.resid))
    #convergence = fit.converged

    #estimate_df = DataFrame(LDR=estimates[1], HDR=estimates[2], ic50=estimates[3],
    #                        slope=estimates[4], rmse=rmse, convergence=convergence);

    #return estimate_df
end

n = 12
dose_range = Array(range(-3, stop=3, length=n));
efficiency_metrics = [100.0, 0.0, 0., 1.0];

R = 10  # Number of repetitions for each condition

rslt_df = DataFrame(sigma=Float64[], repetition=Int[], subset_size=Int[], method=String[], time=Float64[])

synthetic_fn = data_prefix*"/synthetic_data_gaussian_noise_time_$(string(R))_repeats_oni_test.csv";
if isfile(synthetic_fn)
    rm(synthetic_fn)
end


for σ in [1., 10.]

    for r in 1:R
        println("σ=$(string(σ)), repetition $r/$R")
        synthetic_df = make_dataset_with_gaussian_noise(dose_range, efficiency_metrics, σ);
        xs_all = synthetic_df.x
        ys_noisy_all = synthetic_df.y_noisy

        for i in n:-1:2
            xs = @view xs_all[1:i]
            ys_noisy = @view ys_noisy_all[1:i]

            posterior_time = do_inference_with_BiDRA(xs, ys_noisy);
            estimate_time = do_curve_fit_with_LsqFit(xs, ys_noisy);

            push!(rslt_df, (σ, r, i, "BiDRA", posterior_time))
            push!(rslt_df, (σ, r, i, "LsqFit", estimate_time))
        end
    end
end

CSV.write(synthetic_fn, rslt_df)

