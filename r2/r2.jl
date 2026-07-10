using DataFrames, CSV, HDF5, Statistics, StatsBase, Random
using CairoMakie, AlgebraOfGraphics, Colors

include("utils.jl")
include("data_prep.jl")

# dts = DataFrame([
#     (name="Gray", mle="gray", id=:exp_id),
#     (name="gCSI", mle="gCSI", id=:expid), 
#     (name="CTRPv2", mle="ctrpv2", id=:experimentIds)
# ])
# dts.sym = Symbol.(dts.name)

# environ 10 minutes
@time dict_all_df, q_df = prep_all_data(dts, 15, false)
dict_all_dfr, q_dfr = prep_all_data(dts, 15, true)


# Figure

fig_p = prep_figure(dts, dict_all_df, "Pearson", q_df)
save("tmp/figure_3_r2.pdf", fig_p)

fig_s = prep_figure(dts, dict_all_df, "Spearman", q_df)
save("tmp/figure_S1A_r2.pdf", fig_s)

fig_pr = prep_figure(dts, dict_all_dfr, "Pearson", q_dfr)
save("tmp/figure_5_r2.pdf", fig_p)

fig_sr = prep_figure(dts, dict_all_dfr, "Spearman", q_dfr)
save("tmp/figure_S6_r2.pdf", fig_s)

