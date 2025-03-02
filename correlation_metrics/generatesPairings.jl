using DataFrames, HDF5, JLD2

include("../utils.jl")

datasets = ["gray"]#["gCSI", "ctrpv2", "gray"]
info_prefix = "public_datasets/curves_info/"
output_prefix = "public_datasets/"

### Helpful functions
### Remove characters in exp_id name
function replaceChar(df::DataFrame, col::Symbol)
    df[!, col] .= replace.(df[:, col], " " => "_")
    df[!, col] .= replace.(df[:, col], ":" => "_")
    df[!, col] .= replace.(df[:, col], "/" => "_")
    return df
end

function replicatesPairing(df::DataFrame, dt::String) 
    #Remove exact same replicates in gCSI
    if dt == "gCSI"
        idx_replicates = findall(nonunique(df[:, [:Concentration, :Viability, :pairs]]))
        pairs_to_rem = df[idx_replicates, :pairs]
        df = filter(:pairs => x -> x ∉ pairs_to_rem, df)
    end
    
    exp_id = unique(df.exp_id)
    tmp = DataFrame(rep_1=exp_id[1:2:end], rep_2=exp_id[2:2:end])
    
    println("---> Writing pairing h5")
    @time h5open(joinpath(output_prefix, "rep2_pairing.h5"), "r+") do file
        g = create_group(file, dt)
        g["rep_1"] = Array(tmp.rep_1)
        g["rep_2"] = Array(tmp.rep_2)
    end

    return df, tmp
end

function duplicatesAssociation(df::DataFrame, dt::String)
    println("---> Writing list of exp. for R>2")
    @time jldopen(joinpath(output_prefix, "rep_more2_pairing.jld2"), "a+") do file
        file[dt] = combine(groupby(df, :pairs), :exp_id => unique)
    end
end

### Init h5 file
@time h5open(joinpath(output_prefix, "rep2_pairing.h5"), "w") do file
end

@time jldopen(joinpath(output_prefix, "rep_more2_pairing.jld2"), "w") do file
end

### Loop through all 3 datasets
for i in 1:length(datasets)
    dt = datasets[i]
    
    ## Get viabilities and info 
    data_df = getRawData_h5(dt, false)
    info_df = readCSV(info_prefix*dt*"_info.csv", true)[:, 1:3]
    print(names(info_df))
    info_df = replaceChar(info_df, :exp_id)

    ## Create meta df and print info
    meta_df = innerjoin(data_df, info_df, on=[:exp_id => :exp_id])
    meta_df[!, "pairs"] = string.(meta_df.drugid, ":", meta_df.cellid)
    println("Dataset info")
    println(dt, ": ", length(unique(meta_df.cellid)), " cell lines, ", length(unique(meta_df.drugid)), " drugs, ", length(unique(meta_df.pairs)), " pairs")

    ## Count nb x an exp. is replicated
    replicates_df = combine(groupby(meta_df, :pairs), :exp_id=>(length∘unique)=>:n_rep)
    println("Intra-Dataset replicates info")
    println(dt, ": ", sum(replicates_df.n_rep .== 1), " single exp., ", 
                    sum(replicates_df.n_rep .== 2), " duplicates, ",
                    sum(replicates_df.n_rep .> 2), " replicates")

    ## keep replicates with r = 2
    duplicates_list = filter(:n_rep => x -> x == 2, replicates_df)[:, :pairs]
    tmp = sort(filter(:pairs => x -> x ∈ duplicates_list, meta_df), [:pairs, :exp_id, :drugid])
    dup, prg = replicatesPairing(tmp, dt)

    println("Unique pairings info")
    println(dt, ": ", length(unique(dup.pairs)), " R=2 ")

    ## keep replicates with r > 2
    replicates_list = filter(:n_rep => x -> x > 2, replicates_df)[:, :pairs]
    tmp = sort(filter(:pairs => x -> x ∈ replicates_list, meta_df), [:pairs, :exp_id, :drugid])
    duplicatesAssociation(tmp, dt)

    println("\n\n\n")
end