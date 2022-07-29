using DataFrames, HDF5
using CSV
using Optim, GLM, LsqFit
using Glob


function checkFile(fn::String)
    if isfile(fn)
        return fn
    elseif isfile("../$fn")
        return "../$fn"
    else
        return missing
    end
end

function readCSV(fn::String, h::Bool, addExpId::Bool, expId::String) 
    csv_file = DataFrame(CSV.File("$fn", header=h, ntasks=8))
    if addExpId
        csv_file[!,"exp_id"] = repeat([expId], nrow(csv_file))
    end
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

function getRawData_h5(dt::String, localVar::Bool)
    ### Define path
    if localVar
        fn_h5 = checkFile("data/local_$dt"*"_complete.h5")
    else
        fn_h5 = checkFile("data/$dt"*"_complete.h5")
    end

    ### Get list of expID
    expId_list = getExpId_h5(dt)
    ## create StrIndex

    ### Import and merge all responses and concentration
    file = h5open(fn_h5, "r")
    @time expSize = map(e -> size(file[e*"/data"])[1], expId_list)
    

    data = map(e -> read(file, e)["data"], expId_list)
    occ_expId = map(e -> repeat([e], size(read(file, e)["data"])[1]), expId_list)
    data_colName = read(file, "info")["data_colNames"]
    close(file)

    ### Build dataframe
    data_mtx = vcat(data...)
    data_df = DataFrame(data_mtx, :auto)
    rename!(data_df, data_colName)
    data_df[!, :exp_id] = vcat(occ_expId...)
    return data_df
end

function getPosterior_h5(dt::String, localVar::Bool)
    ### Define path
    if localVar 
        fn_h5 = "data/local_$dt"*"_complete.h5"
    else
        fn_h5 = "data/$dt"*"_complete.h5"
    end

    ### Get list of expID
    expId_list = getExpId_h5(dt)

    ### Import and merge all posterior
    file = h5open(fn_h5, "r")
    chains = map(e -> read(file, e)["chains"], expId_list)
    chains_colName = read(file, "info")["chains_colNames"]
    close(file)

    chains_mtx = vcat(chains...)
    chains_df = DataFrame(chains_mtx, :auto)
    rename!(chains_df, chains_colName)

    ### Get experiments listing
    chains_df[!, :exp_id] = repeat(expId_list, inner=4000)
    return chains_df
end

function getPosterior_h5(dt::String, localVar::Bool, expId_list::Array)
    ### Define path
    if localVar 
        fn_h5 = "data/local_$dt"*"_complete.h5"
    else
        fn_h5 = "data/$dt"*"_complete.h5"
    end

    ### Import and merge all posterior
    file = h5open(fn_h5, "r")
    chains = map(e -> read(file, e)["chains"], expId_list)
    chains_colName = read(file, "info")["chains_colNames"]
    close(file)

    chains_mtx = vcat(chains...)
    chains_df = DataFrame(chains_mtx, :auto)
    rename!(chains_df, chains_colName)

    ### Get experiments listing
    chains_df[!, :exp_id] = repeat(expId_list, inner=4000)
    return chains_df
end

function llogistic(param::Array) 
    LDR, HDR, ic50, slope = param    
    return x -> HDR + ((LDR - HDR) / (1 + 10^(slope * (x - ic50))))
end

function getPairings_h5(dt::String)
    fn = checkFile("correlation_metrics/rep2_pairing.h5")
    df = DataFrame(h5read(fn, dt))
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

function getMLestimates(dt::Array,)
    mle_prefix = "data/all_julia_curveFit.csv"
    mle_data = readCSV(mle_prefix, true, false, "")
    mle_data = filter(:dataset => x -> x ∈ dt, mle_data)
    return mle_data
end

function getMLestimates(dt::Array, pairing_df::DataFrame)
    mle_data = getMLestimates(dt)
    
    mle_data = filter(:exp_id => x -> x ∈ pairing_df.rep_1 || x ∈ pairing_df.rep_2, mle_data)
    mle_data = mle_data[:, [:exp_id, :LDR, :HDR, :ic50, :slope, :aac, :dataset, :convergence]]

    mle_tmp = innerjoin(pairing_df, mle_data, on=:rep_1 => :exp_id, renamecols=("" => "_rep1"))
    mle_tmp = innerjoin(mle_tmp, mle_data, on=:rep_2 => :exp_id, renamecols=("" => "_rep2"))
    
    mle_data = filter(row -> !isnan(row.LDR_rep1) && !isnan(row.LDR_rep2), mle_tmp)
    mle_data[!, :aac_rep_1] = replace(mle_data.aac_rep1, Inf => NaN, -Inf => NaN)
    mle_data[!, :aac_rep_2] = replace(mle_data.aac_rep2, Inf => NaN, -Inf => NaN)
    
    return mle_data
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
        f = llogistic(tmp)
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