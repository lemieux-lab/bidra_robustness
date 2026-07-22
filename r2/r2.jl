using DataFrames, CSV, HDF5, JLD2, Statistics, StatsBase, Random
using CairoMakie, AlgebraOfGraphics, Colors

include("utils.jl")
include("data_prep.jl")

cache_fn = "public_datasets/cache_all.jld2"
if isfile(cache_fn)
    @load cache_fn dict_all_df q_df dict_all_dfr q_dfr
else
    # about 12 minutes on 120 cores
    dict_all_df, q_df = prep_all_data(dts, 15, false)
    dict_all_dfr, q_dfr = prep_all_data(dts, 15, true)

    @save cache_fn dict_all_df q_df dict_all_dfr q_dfr
end


# Figure

fig_p = prep_figure(dts, dict_all_df, "Pearson", q_df)
save("tmp/figure_3_r2.pdf", fig_p)

fig_s = prep_figure(dts, dict_all_df, "Spearman", q_df)
save("tmp/figure_S1A_r2.pdf", fig_s)

fig_pr = prep_figure(dts, dict_all_dfr, "Pearson", q_dfr)
save("tmp/figure_5_r2.pdf", fig_p)

fig_sr = prep_figure(dts, dict_all_dfr, "Spearman", q_dfr)
save("tmp/figure_S6_r2.pdf", fig_s)

