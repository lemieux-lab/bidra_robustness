using DataFrames
using Statistics, StatsBase
using ProgressBars

include("../utils.jl")
fn = "_generated_data/viabCorrelations.csv"

### Define dataset to analyze
dt = ARGS[1]

### Get data
println("1. Get data")
@time expId_list = getExpId_h5(dt);
@time si = StrIndex(expId_list);
pairings_df = getPairings_h5(dt, si)
pairings_df.pairID = collect(1:nrow(pairings_df))
viability_df = getRawData_h5(dt, false, si)


println("2. Pair experiments responses")
rep1 = innerjoin(pairings_df[:, [:rep_1, :pairID]], viability_df, on=:rep_1=>:exp_id)
rep2 = innerjoin(pairings_df[:, [:rep_2, :pairID]], viability_df, on=:rep_2=>:exp_id)
pair_shared_dose = innerjoin(rep1, rep2, on=[:pairID, :Concentration], makeunique=true)

### Correlation Analysis
function doCorrelation(df::DataFrame, colNames::Array,  dt::String, description::String, addHeader::Bool)
    r1, r2 = colNames
    corr_df = correlationAnalysis(df[:,r1], df[:,r2])
    corr_df[:, :N] = [length(unique(df.rep_1))]
    corr_df[:, :n] = [nrow(df)]
    corr_df[:, :dataset] = [dt]
    corr_df[:, :description] = [description]
    
    if addHeader 
        CSV.write(fn, corr_df, delim=",", append=false, header=["slope","intercept","r²","rₛ","r","N","n","dataset","description"])
    else
        CSV.write(fn, corr_df, delim=",", append=true)
    end
end

println("3. Do correlation analysis")
println("---> all")
## Considering all pairings
doCorrelation(pair_shared_dose, [:Viability, :Viability_1], dt, "all possible pairing", true)

## Mean response
tmp_mean = combine(groupby(pair_shared_dose, [:rep_1, :Concentration]), :Viability => mean => :mean_1, :Viability_1 => mean => :mean_2)

if nrow(tmp_mean) != nrow(pair_shared_dose)
    println("---> mean")
    doCorrelation(tmp_mean, [:mean_1, :mean_2], dt, "mean per dose", false)

    println("---> bootstrap")
    ## Boostrap
    N = length(unique(tmp_mean.rep_1))
    n = nrow(tmp_mean)
    R = 10000
    rep_corr_df = DataFrame()

    for r in ProgressBar(1:R)
        tmp = mapreduce(g -> DataFrame(g[rand(1:nrow(g)), :]), vcat, groupby(pair_shared_dose, [:rep_1, :Concentration]))
        corr_df = correlationAnalysis(tmp[:,:Viability], tmp[:,:Viability_1])
        append!(rep_corr_df, corr_df)
    end

    function saveRepStats(df::DataFrame, s::Symbol)
        tmp = DataFrame(slope=[], intercept=[] ,r²=[], rₛ=[], r=[], N=[], n=[], dataset=[], description=[])
        push!(tmp, vcat(df[!, s], [N, n, dt, "$R rep $s"]))
        CSV.write(fn, tmp, delim=",", append=true)
    end

    stats_rep = describe(rep_corr_df, :mean, :median, std => :std)
    saveRepStats(stats_rep, :mean)
    saveRepStats(stats_rep, :median)
    saveRepStats(stats_rep, :std)
else 
    println("There are no dose-replicate")
end
