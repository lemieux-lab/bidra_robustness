using DataFrames

include("../utils.jl")

datasets = ["gCSI", "ctrpv2", "gray"]
info_prefix = "/home/golem/scratch/labellec/_DATA/curves_info/"
data_prefix = "/home/golem/scratch/labellec/_DATA/"
output_prefix = "./"

### Dataset-specific
info_expLabels = [:expid, :experimentIds, :exp_id]

### Helpful functions
### Remove characters in exp_id name
function replaceChar(df, col)
    df[!, col] .= replace.(df[:, col], " " => "_")
    df[!, col] .= replace.(df[:, col], ":" => "_")
    df[!, col] .= replace.(df[:, col], "/" => "_")
    return df
end

function replicatesPairing(df, dt) 
    #Remove exact same replicates in gCSI
    if dt == "gCSI"
        idx_replicates = findall(nonunique(df[:, [:Concentration, :Viability, :pairs]]))
        pairs_to_rem = df[idx_replicates, :pairs]
        df = filter(:pairs => x -> x ∉ pairs_to_rem, df)
    end
    
    exp_id = unique(df.exp_id)
    tmp = DataFrame(rep_1=exp_id[1:2:end], rep_2=exp_id[2:2:end])
    
    println("---> Writing pairing files")
    CSV.write(output_prefix*dt*"_rep2_pairing.csv", tmp)
    
    return df, tmp
end

### Loop through all 3 datasets
for i in 1:length(datasets)
    dt = datasets[i]

    ## Get viabilities and info 
    data_df = getRawData([dt], data_prefix, true)
    #data_df = readCSV(data_prefix*dt*"_selected_curves_all.csv", true)
    #data_df[!, :Concentration] = log10.(data_df[:, :Concentration])
    info_df = readCSV(info_prefix*dt*"_info.csv", true)[:, 1:3]
    info_df = replaceChar(info_df, info_expLabels[i])

    ## Create meta df and print info
    meta_df = innerjoin(data_df, info_df, on=[:exp_id => info_expLabels[i]])
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
    println("\n\n\n")
end