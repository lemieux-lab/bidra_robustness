using DataFrames
using CSV
using Optim, GLM, LsqFit

function readCSV(fn, h) 
    csv_file = CSV.File("$fn", header=h)
    return DataFrame(csv_file); 
end

function replaceChar(df, col)
    df[!, col] .= replace.(df[:, col], " " => "_")
    df[!, col] .= replace.(df[:, col], ":" => "_")
    df[!, col] .= replace.(df[:, col], "/" => "_")
    return df
end

function getRawData(dt, dt_pf, h)
    data_df = DataFrame()

    for i in dt 
        tmp_df = readCSV(dt_pf*i*"_selected_curves_all.csv", h)
        tmp_df[!, :dataset] = repeat([i], nrow(tmp_df))
        tmp_df[!, :Concentration] = log10.(tmp_df[:, :Concentration])
        
        maxThresh = 200
        minThresh = -50
        id_extreme = unique(filter(:Viability => x -> x > maxThresh || x < minThresh, tmp_df)[:, :exp_id])
        tmp_df = filter(:exp_id => x -> x ∉ id_extreme, tmp_df)
        
        append!(data_df, replaceChar(tmp_df, :exp_id))
    end

    ### Print summary
    for i in dt
        tmp = filter(:dataset => x -> x == i, data_df)
        println(i, " : ", length(unique(tmp.exp_id)), " experiments")
    end
    println("------------------------------------")
    println("TOTAL : ", length(unique(data_df.exp_id)), " experiments")
    println()

    return data_df
end

function llogistic(param) 
    LDR, HDR, ic50, slope = param    
    return x -> HDR + ((LDR - HDR) / (1 + 10^(slope * (x - ic50))))
end

function getPairings(dt)
    prefix = "/u/labellec/Desktop/bayesian_dose_response/articles/"
    df = readCSV(prefix*"_DATA/"*dt*"_rep2_pairing.csv", true)
    return df
end

function getMLestimates(dt, paired, pairing_df)
    mle_prefix = "/home/golem/scratch/labellec/_RESULTS/MLE_ESTIMATES/all_julia_curveFit.csv"
    mle_data = readCSV(mle_prefix, true)

    mle_data = filter(:dataset => x -> x ∈ dt, mle_data)
    
    if paired
        mle_data = filter(:exp_id => x -> x ∈ pairing_df.rep_1 || x ∈ pairing_df.rep_2, mle_data)
        mle_data = mle_data[:, [:exp_id, :LDR, :HDR, :ic50, :slope, :aac, :dataset, :convergence]]

        mle_tmp = innerjoin(pairing_df, mle_data, on=:rep_1 => :exp_id, renamecols=("" => "_rep_1"))
        mle_tmp = innerjoin(mle_tmp, mle_data, on=:rep_2 => :exp_id, renamecols=("" => "_rep_2"))
    
        mle_data = filter(row -> !isnan(row.LDR_rep_1) && !isnan(row.LDR_rep_2), mle_tmp)

        ### Print summary
        for i in dt
            tmp = filter(:dataset_rep_1 => x -> x == i, mle_data)
            println(i, " : ", length(unique(tmp.rep_1)), " experiments duplicates")
        end
        println("------------------------------------")
        println("TOTAL : ", length(unique(mle_data.rep_1)), " experiments experiments duplicates")
        println()
    end

    return mle_data
end

function getpGXestimates(dt, paired, pairing_df)
    pgx_prefix = "/home/golem/scratch/labellec/_RESULTS/MLE_PHARMACOGX/MLE_PGx/" * dt * "_estimates.csv"
    pgx_data = readCSV(pgx_prefix, true)
    
    pgx_data = replaceChar(pgx_tmp, :exp_id)
    if paired
        pgx_data = filter(:exp_id => x -> x ∈ pairing_df.rep_1 || x ∈ pairing_df.rep_2, pgx_tmp)
    end

    pgx_data[!, :EC50] = log10.(pgx_data[:, :EC50])
    pgx_data[!, :aac_recomputed] = 100 .* pgx_data[:, :aac_recomputed]
    
    return pgx_data
end

function getBIDRAdiagnotics(dt, paired, pairing_df)
    bidra_prefix = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/"*dt*"_diagnostics.csv"
    bidra_data = readCSV(bidra_prefix, true)
    
    if paired
        bidra_data = bidra_data[:, [:HDR, :LDR, :ic50, :slope, :σ]]
        bidra_data = filter(:exp_id => x -> x ∈ pairing_df.rep_1 || x ∈ pairing_df.rep_2, bidra_data)
    end

    return bidra_data
end

function getBIDRAposterior(dt, expId_list)
    posterior_prefix = "/home/golem/scratch/labellec/_RESULTS/"*dt*"_julia_process_all/"
    
    function getPosterior(exp) 
        return readCSV(posterior_prefix*exp*".csv", true)
    end
    
    data_posterior = DataFrame(HDR=[], LDR=[], ic50=[], slope=[], σ=[], exp_id=[])
    
    for e in expId_list
        tmp = getPosterior(e)[:, [:HDR, :LDR, :ic50, :slope, :σ]]
        tmp[!,"exp_id"] = repeat([e], nrow(tmp))
        append!(data_posterior, tmp)
    end
    
    return data_posterior
    
end

function getPairedPosterior(dt, pairing_df)
    function getPosterior(pair, dt)
        results_path = "/home/golem/scratch/labellec/_RESULTS/"*dt*"_julia_process_all/"
        exp1, exp2 = pair
    
        posterior1 = readCSV(results_path*exp1*".csv", true)[:, [:HDR, :LDR, :ic50, :slope, :σ]]
        posterior2 = readCSV(results_path*exp2*".csv", true)[:, [:HDR, :LDR, :ic50, :slope, :σ]]
    
        posterior1[!,:exp_id] = repeat([exp1], nrow(posterior1))
        posterior2[!,:exp_id] = repeat([exp2], nrow(posterior2))
    
        return posterior1, posterior2
    end

    df_posterior = DataFrame(param=[], exp1_val=[], exp2_val=[], exp_id1=[], exp_id2=[])
    
    for i in 1:nrow(pairing_df)
        posterior1, posterior2 = getPosterior(pairing_df[i,:], dt)

        tmp_sampled1 = stack(posterior1[1:4000,:], Not("exp_id"))
        tmp_sampled2 = stack(posterior2[1:4000,:], Not("exp_id"))

        posterior_sampled = DataFrame(param=tmp_sampled1.variable, exp1_val=tmp_sampled1.value, 
            exp2_val=tmp_sampled2.value, exp_id1=tmp_sampled1.exp_id, exp_id2=tmp_sampled2.exp_id)
        append!(df_posterior, posterior_sampled)
    end
    
    return df_posterior 
end

function correlationAnalysis(X, Y)
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

    return DataFrame(slope=lr_slope, r²=rSquared, rₛ=spearman, r=pearson)
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