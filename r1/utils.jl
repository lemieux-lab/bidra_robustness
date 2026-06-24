using DataFrames, HDF5, JLD2, CSV

h5safe(str) = replace(str, ' ' => '_', ':' => '_', '/' => '_')

function identify_replicates(h, dt, id_col)
    k = Set(keys(h) .|> Symbol)
    info_df = CSV.read("public_datasets/curves_info/$(dt)_info.csv", DataFrame; pool = true)
    info_df.cond = [Symbol("$(row.cellid):$(row.drugid)") for row in eachrow(info_df)]
    gdf = groupby(info_df, :cond)
    paired_cond = subset(combine(gdf, nrow => :n), :n => v -> v .== 2)

    rep_1 = Symbol[]
    rep_2 = Symbol[]
    for (i, row) in eachrow(paired_cond) |> enumerate
        ids = gdf[(row.cond,)][!, id_col] .|> h5safe .|> Symbol # assumes ids has length 2
        if all(p ∈ k for p in ids)
            push!(rep_1, ids[1])
            push!(rep_2, ids[2])
        end
    end
    df = DataFrame(rep_1 = rep_1, rep_2 = rep_2)

    return df
end

# This is an approximation of pearson when n becomes infinite
function r_swap(a, b)
    n = length(a)

    sa = 0.0
    sb = 0.0
    saa_sbb = 0.0
    sab = 0.0

    for i in eachindex(a, b)
        x = a[i]
        y = b[i]

        sa += x
        sb += y
        saa_sbb += x*x + y*y
        sab += x*y
    end

    μ = (sa + sb) / (2n)

    num = sab / n - μ^2
    den = saa_sbb / (2n) - μ^2

    return num / den
end


function r_swap_mc(a, b; B=10_000)
    n = length(a)

    μz = (sum(a) + sum(b)) / (2n) # average mid-point
    
    Vz = 0.0
    D  = 0.0
    
    d  = similar(a)
    zd = similar(a)

    for i in eachindex(a, b)
        zi = (a[i] + b[i]) / 2
        zc = (zi - μz)

        d[i] = (a[i] - b[i]) / 2
        zd[i] = zc * d[i]

        Vz += zc * zc
        D  += d[i] * d[i]
    end

    Vz /= n
    D  /= n

    out = 0.0

    for k in 1:B
        ϵ = rand([-1., +1.], n)

        η  = sum(ϵ .* d) / n
        C = sum(ϵ .* zd) / n


        out += (Vz - D + η^2) / sqrt((Vz + D - η^2)^2 - 4C^2)
    end

    return out / B
end

