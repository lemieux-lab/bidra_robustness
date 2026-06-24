using DataFrames, CSV, HDF5, Statistics, StatsBase
using CairoMakie, AlgebraOfGraphics

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

subs = Dict(
    :both_c => all(is_complete; dims=2) |> vec,
    :both_i => all(.~is_complete; dims=2) |> vec,
    :all => trues(nrow(df))
)

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

# Launch analysis

metrics = [:LDR, :HDR, :ic50, :slope]
sub_labels = [:all, :both_c, :both_i]

gen_qs(k) = ((1:k) .- 0.5) ./ k
qs = Dict(Symbol(q) => q for q in gen_qs(5))
tmp_q = gen_qs(7)
q_df = DataFrame(sym=Symbol.(tmp_q), q=tmp_q) #, color=)
push!(q_df, (sym=:med, q=0.5))

res = Dict{Tuple{Symbol, Symbol}, Vector{Float32}}()
res_lock = ReentrantLock()

Threads.@threads for metric ∈ metrics
    for sub_label ∈ sub_labels
        tmp = do_all(mi[subs[sub_label],:,:], mj[subs[sub_label],:,:], chains_colName[metric])
        lock(res_lock) do 
            res[(metric, sub_label)] = tmp
        end
    end
end

all_df = DataFrame()

for metric ∈ metrics, sub_label ∈ sub_labels
    for row ∈ eachrow(q_df)
        row.sym == :med && continue
        push!(all_df, (method=:bidra_quantile, metric=metric, sub=sub_label, q=row.sym, r_swap=quantile(res[(metric, sub_label)], row.q)))
    end

    push!(all_df, (method=:mle, metric=metric, sub=sub_label, q=:med, r_swap=r_swap_mc(di[subs[sub_label], metric], dj[subs[sub_label], metric], 10_000)))
end

# Figure

using AlgebraOfGraphics, ColorSchemes

df_mle = subset(all_df, :method => ByRow(==(:mle)))
df_bidra = subset(all_df, :method => ByRow(!=(:mle)))

method_map = :method => renamer([:mle => "L-M", :bidra_quantile => "BiDRA q", :bidra_hdi => "BiDRA HDI",]) => "Methods"

sub_map = :sub => renamer([:all => "All pairs", :both_c => "Complete pairs\nSD ≥ 20", :both_i => "Incomplete pairs\nSD < 20"])

base_mapping = mapping(
    method_map,
    :r_swap => "Swap-invariant coefficient of concordance",
    row = :metric,
    col = sub_map,
    color = :q,
)

spec_bar_bidra = data(df_bidra) * base_mapping *
    mapping(dodge = :q => sorter(q_df.sym)) *
    visual(BarPlot; dodge_gap = 0.01)

spec_bar_mle = data(df_mle) * base_mapping * visual(BarPlot; gap=1.5)

spec_hl = data((pos = [0.5, 0.75],)) * mapping(:pos) * visual(HLines; linestyle = :dash)

edge_weight(x, γ=1; mid=0.0, edge=1.0) = mid + (edge - mid) * (2abs(x - 0.5))^γ

function my_col(x)
    s = edge_weight(x, 0.5; mid=0., edge=1.0)
    h = x < 0.5 ? 0. : 240.
    l = edge_weight(x, 1; mid=0.0, edge=0.9)
    return HSL(h, s, l)
end

q_df.color = my_col.(q_df.q)

draw(
    spec_hl + spec_bar_mle + spec_bar_bidra,
    scales(Color = (; palette = [row.sym => row.color for row ∈ eachrow(q_df)]));
    axis = (;
        width = 200, height = 200,
        xticklabelrotation = pi / 4,
        xgridvisible = false,
        yminorticks = -0.1:0.1:1.0,
        yminorticksvisible = true,
        yminorgridvisible = true,
        yminorgridcolor = "#aaaaaa",
        yminorgridwidth = 0.5,
        limits = (nothing, nothing, -0.1, 1),
    ),
)

