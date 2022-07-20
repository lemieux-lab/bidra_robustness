### BiDRA code will be modified to export posterior results in HDF5 format rather than CSV
### The following code is a temporary quick n easy fix
using DataFrames
using HDF5

include("utils.jl")

dt = ARGS[1]
#dt = "gray"

function getRawData_csv(dt::String)
    data_fn = "data/$dt"*"_selected_curves_all.csv"
    data_df = readCSV(data_fn, true, false, "")
    data_df[!, :Concentration] = log10.(data_df[:, :Concentration])
        
    maxThresh = 200
    minThresh = -50
    id_extreme = unique(filter(:Viability => x -> x > maxThresh || x < minThresh, data_df)[:, :exp_id])
    data_df = filter(:exp_id => x -> x ∉ id_extreme, data_df)
        
    data_df = replaceChar(data_df, :exp_id)

    ### Print summary
    println("TOTAL : ", length(unique(data_df.exp_id)), " experiments")
    println()
    return data_df
end
    
println("################ $dt ################")
println("1. Define path")
fn_h5 = "data/local_$dt"*"_complete.h5"
data_prefix = "data/$dt"*"_julia_process_all"

println("2. Get exp id")
@time data_df = getRawData_csv(dt)
expId_list = unique(data_df.exp_id)

println("3. Open and write h5 file with individual experiments info")
@time h5open(fn_h5, "w") do file
    tmp_chains = DataFrame()
    tmp_data = DataFrame()

    for e in expId_list
        file_name = joinpath(data_prefix, "$e.csv")
        tmp_chains = CSV.read(file_name, DataFrame, ntasks=12)
        tmp_data = filter(:exp_id => x -> x == e, data_df)[:, [:Concentration, :Viability]]

        g = create_group(file, e)
        g["data"] = Array(tmp_data) 
        g["chains"] = Array(tmp_chains)
    end

    g = create_group(file, "info")
    g["exp_id"] = expId_list
    g["dataset"] = dt
    g["data_colNames"] = names(tmp_data)
    g["chains_colNames"] = names(tmp_chains)
end
println("\n\n")
