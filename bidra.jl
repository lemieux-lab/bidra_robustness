using Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("MCMC_models.jl")
include("utils.jl")

### Batch variables
mod = parse(Int64, ARGS[1])
batch = parse(Int64, ARGS[2])

### Global variables
datasets = ["gCSI"]
data_prefix = "/home/golem/scratch/labellec/_DATA/"
result_prefix = "/home/golem/scratch/labellec/_RESULTS/gCSI_julia_process_all/"
figure_prefix = "/home/golem/scratch/labellec/_RESULTS/gCSI_julia_process_all/FIGURES/"
diagnostic_fn = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/gCSI_diagnostics.csv"
diagnosticTMP_fn = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/TMP_diagnostics"*string(batch)*".csv"
batchTiming_fn = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/batch_timing.csv"

### NUT-s parameters
nChain = 4
nIte = 1000
nAdapt = 1000
δ = 0.65
#Turing.setadbackend(:reversediff)

data_df = getRawData(datasets, data_prefix, true)
expId = unique(data_df.exp_id)

### Create batches and select one
sort!(expId)
subset_expId = [expId[i] for i in 1:length(expId) if i%mod == batch]

##### Subset 
T = @elapsed for e in subset_expId
    println(e)
    subset_df = filter(:exp_id => x -> x == e, data_df)

    ### inference
    println("DO INFERENCE")
    bidra_model = BIDRA(subset_df.Concentration, subset_df.Viability)
    bidra_sampler = NUTS(nAdapt, δ)
    t = @elapsed bidra_chains = sample(bidra_model, bidra_sampler, MCMCThreads(), nIte, nChain)

    ### save posterior values
    println("SAVE POSTERIOR VALUES")
    posterior_df = DataFrame(bidra_chains)
    toSave = posterior_df[:,[:HDR, :LDR, :ic50, :slope, :σ, :chain]]
    CSV.write(result_prefix*e*".csv", toSave)

    ### add to diagnostics
    println("UPDATE DIAGNOTICS")
    medVal = median.(eachcol(toSave[:,[:HDR, :LDR, :ic50, :slope, :σ,]]))
    diagnostics = DataFrame(exp_id=e, batch=batch, time=t, 
                            HDR=medVal[1], LDR=medVal[2], ic50=medVal[3], slope=medVal[4],
                            σ=medVal[5])
    CSV.write(diagnosticTMP_fn, diagnostics, delim=",", append=true)

    ### generated figures
    println("GENERATE FIGURE")
    xDose = minimum(subset_df.Concentration)-1:0.1:maximum(subset_df.Concentration)+1
    X = zeros(nrow(posterior_df),length(xDose))

    curves_df = DataFrame()
    for c in xDose 
        curves_df[!,Symbol(c)] = [0.0]
    end

    for i = 1:nrow(posterior_df)
        tmp = posterior_df[i,[:LDR, :HDR, :ic50, :slope]]
        f = llogistic(tmp)
        y_tmp = f.(xDose)
        
        push!(curves_df, y_tmp)
    end

    med_curve = median.(eachcol(curves_df))
    upper_curve = percentile.(eachcol(curves_df), 97.5)
    lower_curve = percentile.(eachcol(curves_df), 2.5)

    p = Gadfly.plot(layer(x=xDose, y=med_curve, Geom.line()),
                    layer(subset_df, x=:Concentration, y=:Viability, Geom.point()),
                    layer(x=xDose, ymin=lower_curve, ymax=upper_curve, Geom.ribbon()))
    #display(p)
    draw(PDF(figure_prefix*e*".pdf", 4inch, 3inch), p)
end

timing_df = DataFrame(batch=batch, nCurves=length(subset_expId), totalTime=T, dataset=datasets[1])
CSV.write(batchTiming_fn, timing_df, delim=",", append=true)