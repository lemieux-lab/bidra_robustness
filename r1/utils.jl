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


function r_swap_mc(a::AbstractVector{T}, b::AbstractVector{T}, B=10_000) where {T <: AbstractFloat}
    n = length(a)
    half = inv(T(2))
    invn = inv(T(n))

    μz = (sum(a) + sum(b)) * half * invn # average mid-point
    # println(typeof(μz))
    
    Vz = zero(T)
    D = zero(T)
    
    d  = similar(a)
    zd = similar(a)

    for i in eachindex(a, b)
        zi = (a[i] + b[i]) * half
        zc = (zi - μz)

        d[i] = (a[i] - b[i]) * half
        zd[i] = zc * d[i]

        Vz += zc * zc
        D  += d[i] * d[i]
    end

    Vz /= n
    D  /= n

    out = zeros(T, B)

    for k in 1:B
        η  = zero(T)
        C = zero(T)
        # ϵ = rand([one(T), -one(T)], n)

        @simd for i ∈ eachindex(d, zd)
        #     # ϵ = one(T) # rand(Bool) ? -one(T) : one(T)
        #     # ϵ = rand([-one(T), one(T)])
        #     # ϵ = ifelse(rand(Bool), -one(T), one(T))
        #     # ϵ = rand([T(-1), T(+1)])
            if rand(Bool)
                η += d[i]
                C += zd[i]
            else
                η -= d[i]
                C -= zd[i]
            end
        end


        η *= invn
        C *= invn

        out[k] = (Vz - D + η^2) / sqrt((Vz + D - η^2)^2 - 4C^2)
    end

    return mean(out)
end

