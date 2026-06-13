using DataFrames, CSV, HDF5, Statistics, StatsBase

dt = "ctrpv2"#"gCSI"#"gray"
df = DataFrame(h5read("correlation_metrics/rep2_pairing.h5", dt))

# chains_colName = h["info/chains_colNames"][:]

mle_data = CSV.read("public_datasets/all_julia_curveFit.csv", DataFrame; pool = true)

nt(row) = (LDR=row.LDR, HDR=row.HDR, ic50=row.ic50, slope=row.slope, aac=row.aac, rmsd=row.rmsd, convergence=row.convergence)
mle = Dict(row.exp_id => nt(row) for row ∈ eachrow(mle_data) if row.dataset == dt)

# h = h5open("/scratch/lemieuxs/ctrpv2_complete.h5")
# (i, j) = df[1,:]
# # for (i, j) in eachrow(df)
#     println("$i, $j")
#     mi = h["/$i/chains"][:,:]
#     mj = h["/$j/chains"][:,:]
# #     break
# # end

# cor(mi[:,1], mj[:,1])

di = DataFrame(mle[i] for i ∈ df[:,1])
dj = DataFrame(mle[i] for i ∈ df[:,2])

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

cp = acor(di.HDR, dj.HDR, 10000, false)
