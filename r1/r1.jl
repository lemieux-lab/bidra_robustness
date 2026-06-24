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

mi = zeros(nrow(df), l_chain, l_param)
mj = zeros(nrow(df), l_chain, l_param)
buf = zeros(l_chain, l_param)

function copy_and_sort(h, df, buf, m, i, j)
    copyto!(buf, h["/$(df[i,j])/chains"])
    sort!(buf, dims=1)
    m[i, :, :] .= buf
end

@time for (i, row) in enumerate(eachrow(df))
    copy_and_sort(h, df, buf, mi, i, 1)
    copy_and_sort(h, df, buf, mj, i, 2)
end

per_quantile = [r_swap(mi[:, i, chains_colName[:LDR]], mj[:, i, chains_colName[:LDR]]) for i ∈ 1:4000]

a = mi[:, 2000, chains_colName[:LDR]]
b = mj[:, 2000, chains_colName[:LDR]]

res = [r_swap_mc(mi[:, i, chains_colName[:LDR]], mj[:, i, chains_colName[:LDR]]) for i ∈ 1:4000]

i = 2000
@time r_swap_mc(a, b)


r_swap(a, b)
