using DataFrames, CSV, HDF5, Statistics, StatsBase

dt = "ctrpv2"#"gCSI"#"gray"
df = DataFrame(h5read("correlation_metrics/rep2_pairing.h5", dt))
all_id = [df.rep_1; df.rep_2] |> unique

chains_colName = h["info/chains_colNames"][:];

mle_data = CSV.read("public_datasets/all_julia_curveFit.csv", DataFrame; pool = true)

nt(row) = (LDR=row.LDR, HDR=row.HDR, ic50=row.ic50, slope=row.slope, aac=row.aac, rmsd=row.rmsd, convergence=row.convergence)
mle = Dict(row.exp_id => nt(row) for row ∈ eachrow(mle_data) if row.dataset == dt)

h = h5open("public_datasets/bidra/$(dt)_complete.h5")

is_complete = falses(size(df))
for (i, row) in enumerate(eachrow(df))
    is_complete[i, :] .= (std(h["$(row.rep_1)/data"][:,2]) >= 20., std(h["$(row.rep_2)/data"][:,2]) >= 20.)
end

function acor(a, b, n, spearman=false)
    mat = [a b]
    
    if spearman
        mat[:, 1] = tiedrank(mat[:, 1])
        mat[:, 2] = tiedrank(mat[:, 2])
    end
    
    c = zeros(n)
    for i=1:n
        for r in axes(mat, 1)
            if rand(Bool)
                mat[r, 1], mat[r, 2] = mat[r, 2], mat[r, 1]
            end
        end
        @views c[i] = cor(mat[:,1], mat[:,2])
    end
    
    return mean(c), std(c)
end

di = DataFrame(mle[i] for i ∈ df[:,1])
dj = DataFrame(mle[i] for i ∈ df[:,2])

both_c = all(is_complete; dims=2)
both_i = all(.~is_complete; dims=2)

cp = acor(di.HDR[both_c], dj.HDR[both_c], 10000, false)

# To extract the chains...
(i, j) = df[1,:]
l_chain, l_param = size(h["/$(df[1,1])/chains"])

mi = zeros(nrow(df) * l_chain, l_param)
mj = zeros(nrow(df) * l_chain, l_param)
buf = zeros(l_chain, l_param)

for (i, row) in enumerate(eachrow(df))
    println("$i, $(row.rep_1)")
    offset = l_chain * (i-1)
    copyto!(buf, h["/$(df[i,1])/chains"])
    mi[(1:4000) .+ offset, :] .= buf
    copyto!(buf, h["/$(df[i,2])/chains"])
    mj[(1:4000) .+ offset, :] .= buf
end

