using DataFrames, CSV, HDF5, Statistics, StatsBase
include("utils.jl")

dt = "CTRPv2"#"gCSI"#"gray"
h = h5open("public_datasets/bidra/$(dt)_complete.h5")

df = identify_replicates(h, "CTRPv2", :experimentIds)
# identify_replicates(hp, "gCSI", :expid)
# identify_replicates(hp, "Gray", :exp_id)
all_id = [df.rep_1; df.rep_2] |> unique

mle_data = CSV.read("public_datasets/all_julia_curveFit.csv", DataFrame; pool = true)

nt(row) = (LDR=row.LDR, HDR=row.HDR, ic50=row.ic50, slope=row.slope, aac=row.aac, rmsd=row.rmsd, convergence=row.convergence)
mle = Dict(Symbol(row.exp_id) => nt(row) for row ∈ eachrow(mle_data) if row.dataset == lowercase(dt))

chains_colName = Dict(Symbol(s) => i for (i, s) ∈ enumerate(h["info/chains_colNames"][:]))

is_complete = falses(size(df))
for (i, row) in enumerate(eachrow(df))
    is_complete[i, :] .= (std(h["$(row.rep_1)/data"][:,2]) >= 20., std(h["$(row.rep_2)/data"][:,2]) >= 20.)
end

di = DataFrame(mle[i] for i ∈ df[:,1])
dj = DataFrame(mle[i] for i ∈ df[:,2])

both_c = all(is_complete; dims=2) |> vec
both_i = all(.~is_complete; dims=2) |> vec

# To extract the chains...
l_chain, l_param = size(h["/$(df[1,1])/chains"])

mi = zeros(Float32, nrow(df), l_chain, l_param)
mj = zeros(Float32, nrow(df), l_chain, l_param)
buf = zeros(l_chain, l_param)

function copy_and_sort(h, df, buf, m, i, j)
    copyto!(buf, h["/$(df[i,j])/chains"])
    sort!(buf, dims=1)
    m[i, :, :] .= Float32.(buf)
end

@time for (i, row) in enumerate(eachrow(df))
    copy_and_sort(h, df, buf, mi, i, 1)
    copy_and_sort(h, df, buf, mj, i, 2)
end

function do_all(mi, mj, metric_idx)
    res = zeros(eltype(mi), size(mi, 2))
    Threads.@threads for i ∈ axes(mi, 2)
        res[i] = r_swap_mc(mi[:, i, metric_idx], mj[:, i, metric_idx], 10_000)
    end
    return res
end

@time res = do_all(mi, mj, chains_colName[:LDR])
quantile(res, [0.025, 0.25, 0.5, 0.75, 0.975])

