push!(LOAD_PATH, "Utils/")
using DataFrames, HDF5
using CSV
using Optim, GLM, LsqFit
using Glob
using ProgressBars
using Utils ### from lemieux-lab


function checkFile(fn::String)
    if isfile(fn)
        return fn
    elseif isfile("../$fn")
        return "../$fn"
    else
        return missing
    end
end

function readCSV(fn::String, h::Bool) 
    csv_file = DataFrame(CSV.File("$fn", header=h, ntasks=8))
    return csv_file 
end

function readCSV(fn::String, h::Bool, expId::String, si::StrIndex) 
    csv_file = DataFrame(CSV.File("$fn", header=h, ntasks=8))
    csv_file[!,"exp_id"] = repeat([si.str2id[expId]], nrow(csv_file))
    return csv_file 
end

function replaceChar(df::DataFrame, col::Symbol)
    df[!, col] .= replace.(df[:, col], " " => "_")
    df[!, col] .= replace.(df[:, col], ":" => "_")
    df[!, col] .= replace.(df[:, col], "/" => "_")
    return df
end

function getExpId_h5(dt::String)
    fn_h5 = checkFile("data/$dt"*"_complete.h5")
    return h5read(fn_h5, "info")["exp_id"]
end

function getRawData_h5(dt::String, localVar::Bool, si::StrIndex)
    ### Define path
    if localVar
        fn_h5 = checkFile("data/local_$dt"*"_complete.h5")
    else
        fn_h5 = checkFile("data/$dt"*"_complete.h5")
    end

    ### Get list of expID
    expId_list = getExpId_h5(dt)

    ### Import responses and concentration
    file = h5open(fn_h5, "r")
    expSize = map(e -> size(file[e*"/data"])[1], expId_list)
    size_tot = sum(expSize)
    
    ## Alocate memory for each column
    concentration_list = Array{Float32, 1}(undef, size_tot)
    viability_list = Array{Float32, 1}(undef, size_tot)
    id_list = Array{Int32, 1}(undef, size_tot)
    pos = 1

    for e in ProgressBar(expId_list)
        n = expSize[findfirst(x -> x == e, expId_list)]
        tmp = read(file, e)["data"]

        concentration_list[pos:pos+n-1] = tmp[1:n]
        viability_list[pos:pos+n-1] = tmp[n+1:end]
        id_list[pos:pos+n-1] = repeat([si.str2id[e]], n)
        
        pos += n
    end
    close(file)

    ### Build dataframe
    data_df = DataFrame(Concentration=concentration_list, Viability=viability_list, exp_id=id_list, dataset=repeat([dt], length(concentration_list)))
    return data_df
end

function getPosterior_h5(dt::String, localVar::Bool, si::StrIndex)
    ### Define path
    if localVar 
        fn_h5 = "data/local_$dt"*"_complete.h5"
    else
        fn_h5 = "data/$dt"*"_complete.h5"
    end

    ### Import all posterior
    expId_list = getExpId_h5(dt)
    file = h5open(fn_h5, "r")

    ## Alocate memory for each column
    chains_colName = read(file, "info")["chains_colNames"]
    expSize = map(e -> size(file[e*"/chains"])[1], expId_list)
    size_tot = sum(expSize)

    chains_mtx = Array{Float32, 2}(undef, size_tot, length(chains_colName))
    id_list = Array{Int32, 1}(undef, size_tot)
    pos = 1

    for e in ProgressBar(expId_list)
        n = expSize[findfirst(x -> x == e, expId_list)]
        tmp = read(file, e)["chains"]

        chains_mtx[pos:pos+n-1,1:length(chains_colName)]=tmp
        id_list[pos:pos+n-1] = repeat([si.str2id[e]], n)
        
        pos += n
    end
    close(file)

    chains_df = DataFrame(chains_mtx, :auto)
    rename!(chains_df, chains_colName)
    chains_df[!, :exp_id] = id_list
    return chains_df
end

function getPosterior_h5(dt::String, localVar::Bool, si::StrIndex, expId_list::Array)
    ### Define path
    if localVar 
        fn_h5 = "data/local_$dt"*"_complete.h5"
    else
        fn_h5 = "data/$dt"*"_complete.h5"
    end

    ### Import all posterior
    file = h5open(fn_h5, "r")

    ## Alocate memory for each column
    chains_colName = read(file, "info")["chains_colNames"]
    expSize = map(e -> size(file[e*"/chains"])[1], expId_list)
    size_tot = sum(expSize)

    chains_mtx = Array{Float32, 2}(undef, size_tot, length(chains_colName))
    id_list = Array{Int32, 1}(undef, size_tot)
    pos = 1

    for e in ProgressBar(expId_list)
        n = expSize[findfirst(x -> x == e, expId_list)]
        tmp = read(file, e)["chains"]

        chains_mtx[pos:pos+n-1,1:length(chains_colName)]=tmp
        id_list[pos:pos+n-1] = repeat([si.str2id[e]], n)
        
        pos += n
    end
    close(file)

    chains_df = DataFrame(chains_mtx, :auto)
    rename!(chains_df, chains_colName)
    chains_df[!, :exp_id] = id_list
    return chains_df
end

function llogistic(param::Array) 
    LDR, HDR, ic50, slope = param    
    return x -> HDR + ((LDR - HDR) / (1 + 10^(slope * (x - ic50))))
end

function getPairings_h5(dt::String, si::StrIndex)
    fn = checkFile("correlation_metrics/rep2_pairing.h5")
    df = DataFrame(h5read(fn, dt))

    df[!,:rep_1] = [si.str2id[v] for v in df[!,:rep_1]]
    df[!,:rep_] = [si.str2id[v] for v in df[!,:rep_2]]
    return df
end

function getPairedPosterior_h5(dt::String)
    pairings_df = getPairings_h5(dt)

    posterior_rep1 = getPosterior_h5(dt, false, Array(pairings_df.rep_1))
    rename!(posterior_rep1, map(x -> "$x"*"_rep1", names(posterior_rep1)))

    posterior_rep2 = getPosterior_h5(dt, false, Array(pairings_df.rep_2))
    rename!(posterior_rep2, map(x -> "$x"*"_rep2", names(posterior_rep2)))

    return hcat(posterior_rep1, posterior_rep2)
end

function getMLestimates(dt::String, si::StrIndex)
    mle_prefix = "data/all_julia_curveFit.csv"
    mle_data = readCSV(mle_prefix, true)

    ## Only select estimate for datasets
    mle_data_dt = filter(:dataset => x -> x == dt, mle_data)
    mle_data_dt[!, :exp_id] = [si.str2id[v] for v in mle_data_dt.exp_id]
    return mle_data_dt
end

function getMLestimates(dt::String, si::StrIndex, pairing_df::DataFrame)
    mle_data = getMLestimates(dt,st)
    
    mle_data = filter(:exp_id => x -> x ∈ pairing_df.rep_1 || x ∈ pairing_df.rep_2, mle_data)
    mle_data = mle_data[:, [:exp_id, :LDR, :HDR, :ic50, :slope, :aac, :dataset, :convergence]]

    mle_tmp = innerjoin(pairing_df, mle_data, on=:rep_1 => :exp_id, renamecols=("" => "_rep1"))
    mle_tmp = innerjoin(mle_tmp, mle_data, on=:rep_2 => :exp_id, renamecols=("" => "_rep2"))
    
    mle_data_paired = filter(row -> !isnan(row.LDR_rep1) && !isnan(row.LDR_rep2), mle_tmp)
    mle_data_paired[!, :aac_rep_1] = replace(mle_data.aac_rep1, Inf => NaN, -Inf => NaN)
    mle_data_paired[!, :aac_rep_2] = replace(mle_data.aac_rep2, Inf => NaN, -Inf => NaN)
    
    return mle_data_paired
end

function correlationAnalysis(X::Array, Y::Array)
    ### Linear Regression
    data = DataFrame(x=X, y=Y)
    data[!, :x] = convert.(Float64, data[:, :x])
    data[!, :y] = convert.(Float64, data[:, :y])
    lf = lm(@formula(x ~ y), data)

    ### Linear Regression coefficient
    lr_slope = GLM.coef(lf)[2]
    rSquared = r2(lf)
    
    ### Correlation coefficient
    spearman = corspearman(data[:, :x], data[:, :y])
    pearson = cor(data[:, :x], data[:, :y])

    return DataFrame(slope=lr_slope, intercept=GLM.coef(lf)[1], r²=rSquared, rₛ=spearman, r=pearson)
end

function getPosteriorCurves(posterior_df, xmin, xmax) 
    xDose = xmin-1:0.1:xmax+2

    curves_df = DataFrame()
    for c in xDose 
        curves_df[!,Symbol(c)] = [0.0]
    end

    for i = 1:nrow(posterior_df)
        tmp = posterior_df[i,[:LDR, :HDR, :ic50, :slope]]
        f = llogistic(Array(tmp))
        y_tmp = f.(xDose)
        
        push!(curves_df, y_tmp)
    end 
    
    return curves_df
end

###### Compute AAC ####
## param: LDR, HDR, IC50, slope
function computeAAC(dose, viability, params)
    a = minimum(dose)
    b = maximum(dose)
    
    viability = viability ./ 100
    params[1] = params[1] / 100
    params[2] = params[2] / 100

    if params[2] == 1 
        aac = 0
    elseif params[4] == 0 
        aac = (params[1] - params[2]) / 2 
    else 
        Δ = b - a
        upDiv = 1 + 10 ^ (params[4] * (b - params[3]))
        downDiv = 1 + 10 ^ (params[4] * (a - params[3]))
        aac = ((params[1] - params[2]) / (params[4] * params[1] * Δ)) * log10(upDiv / downDiv)
    end

    params[1] = params[1] * 100
    params[2] = params[2] * 100

    return aac*100
end

function get_HDR_prior(N)
    λₖ = [0.4, 0.5, 0.1]
    hdr_mm = MixtureModel([SkewNormal(0, 10, 1), Uniform(0, 100), SkewNormal(100, 20, -5)], λₖ)
    hdr_prior = rand(hdr_mm, N)
    return hdr_prior
end

function get_LDR_prior(N)
    ldr_prior = rand(Normal(100,10), N)
    return ldr_prior
end

function get_ic50_prior(N)
    ic50_prior = rand(Normal(0,10), N)
    return ic50_prior
end

function get_slope_prior(N)
    slope_prior = rand(LogNormal(0.5, 1), N)
    return slope_prior
end