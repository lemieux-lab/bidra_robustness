using Gadfly, StatsPlots
using Cairo, Fontconfig

include("MCMC_models.jl")
include("utils.jl")

### Global variables
datasets = ["gCSI", "ctrpv2","gray"]
data_prefix = "/home/golem/scratch/labellec/_DATA/"

### NUT-s parameters
nChain = 2
nIte = 2000
nAdapt = 1000
δ = 0.65
Turing.setadbackend(:reversediff)


function getMaxConcViab()
    data_df = getRawData(datasets, data_prefix, true)
    maxConcentration_df = combine(groupby(data_df, :exp_id), 
        :Concentration => maximum => :Concentration, 
        [:Concentration, :Viability] => ((a,b) -> b[argmax(a)])  => :Viability,
        [:Concentration, :dataset] => ((a,b) -> b[argmax(a)])  => :dataset)

    return maxConcentration_df
end

function equalSampling(df, n)
    subset_df = DataFrame(exp_id=[], Concentration=[], Viability=[], dataset=[])

    for i in datasets
        tmp = filter(:dataset => x -> x == i, df)
        rand_idx = rand(1: nrow(tmp), n)
        append!(subset_df, tmp[rand_idx, :])
    end

    return subset_df
end

function runAnalysis(HDR_sts, LDR_sts, n)
    inferedMetrics = []

    if HDR_sts 
        maxConcViab_df = getMaxConcViab()
        dataSampled_df = equalSampling(maxConcViab_df, n[1])

        ##Infer parameters
        normalBeta_model = modelNormalBeta(dataSampled_df.Viability)
        normalBeta_sampler = NUTS(nAdapt, δ)
        normalBeta_chains = sample(normalBeta_model, normalBeta_sampler, MCMCThreads(), nIte, nChain)

        modelMedian = [median(normalBeta_chains[:α]), median(normalBeta_chains[:β]), median(normalBeta_chains[:σ])]
        append!(inferedMetrics, [modelMedian])

        x_posterior = vcat(rand.(Normal.(rand(Beta(modelMedian[1], modelMedian[2]), 40000) .* 100, modelMedian[3]), 1)...)

        p = Gadfly.plot(layer(x=x_posterior, color=["Posterior"], Geom.density(bandwidth=0.1)),
                        layer(x=dataSampled_df.Viability, color=["Sampled Data"], Geom.histogram(density=true)),
                        Guide.title("nα="*string(round(modelMedian[1], digits=3))*
                                " β="*string(round(modelMedian[2], digits=3))*
                                " σ="*string(round(modelMedian[3], digits=3))))

        display(p)
    end

    return inferedMetrics
end

benchmark_df = DataFrame(α=Float64[], β=Float64[], σ=Float64[], μ=Float64[], timeSec=Float64[], R=Int64[], N=Int64[], ParamName=[])
benchmark_fn = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/hdr_sampling_inference.csv"
N = 1000
for i in 1:10
    t = @elapsed metricsVal= runAnalysis(true, false, [N, 0])
    tmp = DataFrame(α=metricsVal[1][1], β=metricsVal[1][2], σ=metricsVal[1][3], μ=missing, timeSec=t, R=i, N=N, ParamName="HDR")
    CSV.write(benchmark_fn, tmp, delim=",", append=true)  
end