using DataFrames, CSV, HDF5, Statistics, StatsBase, Random
using CairoMakie, AlgebraOfGraphics, Colors, Format

include("utils.jl")
include("data_prep.jl")

all_pm = Dict{String, PairedMetrics}()
for dt ∈ eachrow(dts)
    all_pm[dt.name] = PairedMetrics(dt)
end


function prep_aac(aac_df, dt, sub, pm)
    keep = pm.subs[sub]

    append!(aac_df, DataFrame(dt=dt, method="Levenberg-Marquardt", sub=sub, di=pm.di.aac[keep], dj=pm.dj.aac[keep]))

    aac = pm.chains_colName[:aac]
    mid = size(pm.mi, 2) ÷ 2

    append!(aac_df, DataFrame(dt=dt, method="BiDRA median", sub=sub, di=pm.mi[keep, mid, aac], dj=pm.mj[keep, mid, aac]))
end

aac_df = DataFrame(dt=String[], method=String[], sub=Symbol[], di=Float64[], dj=Float64[])
for dt ∈ dts.name
    for sub ∈ [:all, :both_c, :both_i]
        prep_aac(aac_df, dt, sub, all_pm[dt])
    end
end


keep = isfinite.(aac_df.di)
keep = keep .&& isfinite.(aac_df.dj)
nrow(aac_df) - sum(keep)
keep = keep .&& (aac_df.di .≥ 0.)
keep = keep .&& (aac_df.dj .≥ 0.)
nrow(aac_df) - sum(keep)
keep = keep .&& (aac_df.di .< 101.)
keep = keep .&& (aac_df.dj .< 101.)
nrow(aac_df) - sum(keep)

keepat!(aac_df, keep)

sub_map = :sub => renamer([:all => "All pairs", :both_c => "Complete pairs\nSD ≥ 20", :both_i => "Incomplete pairs\nSD < 20"])

set_theme!(Theme(
    fonts = (
        regular = "Noto Sans",
        bold = "Noto Sans Bold"
)))

layer_1 =
    data(filter(r -> r.dt == "gCSI", aac_df)) * 
    mapping(
        :di => "AAC Rep1", :dj => "AAC Rep2",
        row= :method => sorter(["Levenberg-Marquardt", "BiDRA median"]), col=sub_map) * 
    visual(Hexbin; cellsize=3, threshold=1, colormap=["gray95", "gray35", "black"])
fig = draw(layer_1; axis = (;
    width=200,
    height=200
)
)

save("tmp/figure_S4A.pdf", fig)


aac_gdf = groupby(aac_df, [:dt, :method, :sub])
aac_cdf = combine(aac_gdf,
    :di => length => :n,
    [:di, :dj] => ((di, dj) -> r_swap_mc(di, dj, 10_000, false)) => :Pearson,
    [:di, :dj] => ((di, dj) -> r_swap_mc(di, dj, 10_000, true)) => :Spearman,
)
sdf = stack(aac_cdf, [:Pearson, :Spearman]; variable_name=:cor_name, value_name=:r_swap_mc)

sdf.value = ["$(format(r.r_swap_mc, precision=2)) (n=$(format(r.n, commas=true)))" for r ∈ eachrow(sdf)]
sdf.col_name = ["$(r.cor_name)_$(r.dt)" for r ∈ eachrow(sdf)]

CSV.write("tmp/figure_S4B.csv", unstack(sdf, [:sub, :method], :col_name, :value))
