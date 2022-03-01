using DataFrames
using CSV

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
    return x -> HDR + (LDR - HDR) / (1 + exp(slope * (x - ic50)))
end

