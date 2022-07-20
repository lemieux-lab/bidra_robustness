using DataFrames
using Statistics, StatsBase
using ProgressBars

include("../utils.jl")
fn = checkFile("_generated_data/viabCorrelations.csv")

### Define dataset to analyze
dt = ARGS[1]
#dt = "gray"

### Get data
println("1. Get data")
pairings_df = getPairings_h5(dt)
pairings_df.pairID = collect(1:nrow(pairings_df))
viability_df = getRawData_h5(dt, false)

println("2. Pair experiments responses")
rep1 = innerjoin(pairings_df[:, [:rep_1, :pairID]], viability_df, on=:rep_1=>:exp_id)
rep2 = innerjoin(pairings_df[:, [:rep_2, :pairID]], viability_df, on=:rep_2=>:exp_id)
pair_shared_dose = innerjoin(rep1, rep2, on=[:pairID, :Concentration], makeunique=true)

### Correlation Analysis
function doCorrelation(df::DataFrame, colNames::Array,  dt::String, description::String)
    r1, r2 = colNames
    corr_df = correlationAnalysis(df[:,r1], df[:,r2])
    corr_df[:, :N] = [length(unique(df.rep_1))]
    corr_df[:, :n] = [nrow(df)]
    corr_df[:, :dataset] = [dt]
    corr_df[:, :description] = [description]
    CSV.write(fn, corr_df, delim=",", append=true)
end

println("3. Do correlation analysis")
println("---> all")
## Considering all pairings
doCorrelation(pair_shared_dose, [:Viability, :Viability_1], dt, "all possible pairing")

## Mean response
tmp = combine(groupby(pair_shared_dose, [:rep_1, :Concentration]), :Viability => mean => :mean_1, :Viability_1 => mean => :mean_2)

if nrow(tmp) != nrow(pair_shared_dose)
    println("---> mean")
    doCorrelation(tmp, [:mean_1, :mean_2], dt, "mean per dose")

    println("---> bootstrap")
    ## Boostrap
    N = length(unique(tmp.rep_1))
    n = nrow(tmp)
    R = 10000
    rep_corr_df = DataFrame()

    for r in ProgressBar(1:R)
        tmp = mapreduce(g -> DataFrame(g[rand(1:nrow(g)), :]), vcat, groupby(pair_shared_dose, [:rep_1, :Concentration]))
        corr_df = correlationAnalysis(tmp[:,:Viability], tmp[:,:Viability_1])
        append!(rep_corr_df, corr_df)
    end

    function saveRepStats(df::DataFrame, s::Symbol)
        tmp = DataFrame(slope=[], intercept=[] ,r²=[], rₛ=[], r=[], N=[], dataset=[], description=[])
        push!(tmp, vcat(df[!, s], [N, n, dt, "$R rep $s"]))
        CSV.write(fn, tmp, delim=",", append=true)
    end

    stats_rep = describe(rep_corr_df, :mean, :median, std => :std)
    saveRepStats(stats_rep, :mean)
else 
    println("There are no dose-replicate")
end
