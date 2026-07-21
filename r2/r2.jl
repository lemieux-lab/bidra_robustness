using DataFrames, CSV, HDF5, Statistics, StatsBase, Random
using CairoMakie, AlgebraOfGraphics, Colors

include("utils.jl")
include("data_prep.jl")

# about 10 minutes on 72 cores
@time dict_all_df, q_df = prep_all_data(dts, 15, false)
@time dict_all_dfr, q_dfr = prep_all_data(dts, 15, true)


# Figure

fig_p = prep_figure(dts, dict_all_df, "Pearson", q_df)
save("tmp/figure_3_r2.pdf", fig_p)

fig_s = prep_figure(dts, dict_all_df, "Spearman", q_df)
save("tmp/figure_S1A_r2.pdf", fig_s)

fig_pr = prep_figure(dts, dict_all_dfr, "Pearson", q_dfr)
save("tmp/figure_5_r2.pdf", fig_p)

fig_sr = prep_figure(dts, dict_all_dfr, "Spearman", q_dfr)
save("tmp/figure_S6_r2.pdf", fig_s)



pm = PairedMetrics(dts[3,:])
q_df = prep_q_df(15)
q_df.color = my_col.(q_df.q)

using BenchmarkTools
@btime r_swap_mc(pm.mi[pm.subs[:all],1,1], pm.mj[pm.subs[:all],1,1], 10_000, false)
@btime spearman_swap_mc(pm.mi[pm.subs[:all],1,1], pm.mj[pm.subs[:all],1,1], 100)

a, b = pm.mi[pm.subs[:all],1,1], pm.mj[pm.subs[:all],1,1]

function spearman_swap_mc2(
    a::AbstractVector{T},
    b::AbstractVector{T},
    B::Integer
) where {T <: AbstractFloat}
    rng = Random.default_rng()
    pooled = [a; b]
    o = sortperm(pooled)
    g = similar(o)
    g_size = Int[]
    n = length(a)

    ngs = 0
    previous = zero(T)

    for (k, idx) in pairs(o)
        value = pooled[idx]

        if k == 1 || value != previous
            ngs += 1
            push!(g_size, 0)
            previous = value
        end

        g[idx] = ngs
        g_size[ngs] += 1
    end

    ga = view(g, 1:n)
    gb = view(g, n+1:2n)

    swap = BitVector(undef, n)
    nx = zeros(Int, ngs)

    rank_x = Vector{T}(undef, ngs)
    rank_y = similar(rank_x)

    μ = T(n + 1) / T(2)
    ρsum = zero(T)

    for _ in 1:B
        rand!(rng, swap)
        fill!(nx, 0)

        # Count pooled-value gs assigned to x.
        for i in eachindex(a, b)
            nx[swap[i] ? ga[i] : gb[i]] += 1
        end

        # Reconstruct exact midranks on both axes.
        cum_x = 0
        cum_y = 0

        for g in eachindex(g_size)
            ngx = nx[g]
            ngy = g_size[g] - ngx

            rank_x[g] = T(cum_x) + T(ngx + 1) / T(2)
            rank_y[g] = T(cum_y) + T(ngy + 1) / T(2)

            cum_x += ngx
            cum_y += ngy
        end

        sxx = zero(T)
        syy = zero(T)
        sxy = zero(T)

        for i in eachindex(a)
            if swap[i]
                rx = rank_x[ga[i]]
                ry = rank_y[gb[i]]
            else
                rx = rank_x[gb[i]]
                ry = rank_y[ga[i]]
            end

            dx = rx - μ
            dy = ry - μ

            sxx += dx * dx
            syy += dy * dy
            sxy += dx * dy
        end

        denominator = sqrt(sxx * syy)

        ρsum += sxy / denominator
    end

    return ρsum / T(B)

end

@btime spearman_swap_mc2(a, b, 10_000)