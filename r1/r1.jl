using DataFrames, CSV, HDF5, Statistics, StatsBase
using CairoMakie, AlgebraOfGraphics, Colors

include("utils.jl")

dts = DataFrame([
    (name="CTRPv2", mle="ctrpv2", id=:experimentIds), 
    (name="gCSI", mle="gCSI", id=:expid), 
    (name="Gray", mle="gray", id=:exp_id)
])
dts.sym = Symbol.(dts.name)

q_df = prep_q_df(7)

all_df = DataFrame()
for dt ∈ eachrow(dts)
    prep_data!(dt, q_df, all_df)
end


# Figure

function dataset_spec(all_df, dataset; show_y=false, show_row=false)
    df_mle = subset(all_df, :dataset => ByRow(==(dataset)), :method => ByRow(==(:mle)))
    df_bidra = subset(all_df, :dataset => ByRow(==(dataset)), :method => ByRow(!=(:mle)))

    method_map = :method => renamer([:mle => "Levenberg-\nMarquardt", :bidra_quantile => "BiDRA posteriors\nper quantile", :bidra_hdi => "BiDRA HDI",]) => "Methods"
    sub_map = :sub => renamer([:all => "All pairs", :both_c => "Complete pairs\nSD ≥ 20", :both_i => "Incomplete pairs\nSD < 20"])
    row_map = show_row ?
        (:metric => renamer([:LDR => "LDR", :HDR => "HDR", :ic50 => rich"IC$_{50}$", :slope => "Slope"])) :
        (:metric => renamer([:LDR => "", :HDR => "", :ic50 => "", :slope => ""]))

    base_mapping = mapping(
        method_map,
        :r_swap => "Swap-invariant coefficient of concordance",
        row = row_map,
        col = sub_map,
        color = :q,
    )

    spec_bar_bidra = data(df_bidra) * base_mapping *
        mapping(dodge = :q => sorter(q_df.sym)) *
        visual(BarPlot; dodge_gap = 0.01)

    spec_bar_mle = data(df_mle) * base_mapping * visual(BarPlot; gap=1.5)

    spec_hl = data((pos = [0.5, 0.75],)) * mapping(:pos) * visual(HLines; linestyle=:dash)
    spec_hl += data((pos = [0.],)) * mapping(:pos) * visual(HLines; linestyle=:solid)

    return spec_hl + spec_bar_mle + spec_bar_bidra
end

edge_weight(x, γ=1; mid=0.0, edge=1.0) = mid + (edge - mid) * (2abs(x - 0.5))^γ

function my_col(x)
    s = edge_weight(x, 0.5; mid=0., edge=1.0)
    h = x < 0.5 ? 0. : 240.
    l = edge_weight(x, 1; mid=0.0, edge=0.9)
    return HSL(h, s, l)
end

q_df.color = my_col.(q_df.q)

fig = Figure(size=(1800, 800))

specs = [dataset_spec(all_df, dts[i, :sym], show_y = (i==1), show_row = (i==nrow(dts))) for i ∈ 1:nrow(dts)]

for (i, spec) ∈ enumerate(specs)
    axis = (
        width = 200, height = 200,
        xticklabelrotation = pi / 4,
        xticklabelsize = 18.,
        xgridvisible = false,
        xlabel = "",
        yminorgridvisible = true,
        yminorgridcolor = "#aaaaaa",
        yminorgridwidth = 0.5,
        limits = (nothing, nothing, -0.2, 1),
        titlesize = 20
    )
    if i == 1
        axis = (axis...,
            yticks = [0.0, 0.5, 1.0],
            yticklabelsize = 18.,
            yminorticks = -0.1:0.1:1.0,
            yminorticksvisible = true,
            ylabel = "Swap-invariant coefficient of concordance",
            ylabelsize = 24
        )
    else
        axis = (axis...,
            yticks = Float32[],
            yticksvisible = false,
            yminorticksvisible = false,
            ylabelvisible = false
        )
    end

    Label(
        fig[1, i], dts[i, :name];
        fontsize = 24,
        font = :bold,
        tellwidth = false
    )

    draw!(fig[2, i],
        spec,
        scales(Color = (; palette = [row.sym => row.color for row ∈ eachrow(q_df)], legend = false));
        axis = axis
    )
end

colgap!(fig.layout, 30)

# Label(fig[1, 1, Left()], "A", font = :bold, fontsize = 28, halign = :left)
# Label(fig[1, 2, Left()], "B", font = :bold, fontsize = 28, halign = :left)
# Label(fig[1, 3, Left()], "C", font = :bold, fontsize = 28, halign = :left)

resize_to_layout!(fig)

fig
