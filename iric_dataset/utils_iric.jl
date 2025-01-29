using Glob
using Turing, MCMCChains
using Distributions
using FillArrays
using Logging
using CSV, DataFrames

Logging.disable_logging(Logging.Warn)

function getViability()
    input = "IRIC_anonymized.csv"
    return DataFrame(CSV.File(input))
end

function llogistic(param::Array) 
    LDR, HDR, ic50, slope = param    
    return x -> LDR + ((HDR - LDR) / (1 + 10^(slope * (ic50 - x))))
end

@model function BIDRA(xs::Array, ys::Array) 
    λₖ = [0.4, 0.5, 0.1]
    HDR ~ mm = MixtureModel([SkewNormal(100, 10, 1), Uniform(0, 100), SkewNormal(0, 10, -1)], λₖ)
    
    LDR ~ Normal(0,10)
    ic50 ~ Normal(0,10)
    slope ~ LogNormal(0.5, 1)
    
    σ ~ LogNormal(1, 1)

    for i in 1:length(xs)
        f = llogistic([LDR, HDR, ic50, slope])
        ys[i] ~ Normal(f(xs[i]), σ)
    end
end

function run_BiDRA(subsetDf::DataFrame)
    nChain = 4
    nIte = 1000
    nAdapt = 1000
    δ = 0.65
    Turing.setprogress!(false)

    bidra_model = BIDRA(subsetDf.log10, subsetDf.inhib_col)
    bidra_sampler = NUTS(nAdapt, δ)
    bidra_chains = sample(bidra_model, bidra_sampler, MCMCThreads(), nIte, nChain)
    return(DataFrame(bidra_chains))     
end

function getPosterior(posterior_path::String)
    FILES = glob("*.csv", posterior_path)
    global posterior_df = DataFrame(LDR=Float64[], HDR=Float64[], ic50=Float64[], slope=Float64[], compound=String[])
    for f in FILES
        tmp = CSV.read(f, header=true, DataFrame)
        compound = first(split(last(split(f, "/")), "."))
        
        tmp[!, :compound] = repeat([compound], nrow(tmp))
        global posterior_df = vcat(posterior_df, tmp[:, [:LDR, :HDR, :ic50, :slope, :compound]])
    end

    return posterior_df
end

numstring(x, n0) = string(x + 10^n0)[2:end]