using DataFrames, HDF5, JLD2, CSV

h5safe(str) = replace(str, ' ' => '_', ':' => '_', '/' => '_')

function identify_replicates(h, dt, randomize=false)
    k = Set(keys(h) .|> Symbol)
    info_df = CSV.read("public_datasets/curves_info/$(dt.name)_info.csv", DataFrame; pool = true)
    info_df.cond = [Symbol("$(row.cellid):$(row.drugid)") for row in eachrow(info_df)]
    gdf = groupby(info_df, :cond)
    paired_cond = subset(combine(gdf, nrow => :n), :n => v -> v .== 2)

    rep_1 = Symbol[]
    rep_2 = Symbol[]
    for (i, row) in eachrow(paired_cond) |> enumerate
        ids = gdf[(row.cond,)][!, dt.id] .|> h5safe .|> Symbol # assumes ids has length 2
        if all(p ∈ k for p in ids)
            push!(rep_1, ids[1])
            push!(rep_2, ids[2])
        end
    end
    df = DataFrame(rep_1 = rep_1, rep_2 = rep_2)
    randomize && (df.rep_2 = shuffle(df.rep_2))

    return df
end

function r_swap_mc(a::AbstractVector{T}, b::AbstractVector{T}, B=10_000, spearman=false) where {T <: AbstractFloat}
    spearman && return spearman_swap_mc(a, b, B)

    n = length(a)
    half = inv(T(2))
    invn = inv(T(n))

    if spearman
        a = tiedrank(a)
        b = tiedrank(b)
    end

    μz = (sum(a) + sum(b)) * half * invn # average mid-point
    
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

        @simd for i ∈ eachindex(d, zd)
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

function dataset_spec(all_df, dataset, q_df; show_row=false)
    df_mle = subset(all_df, :dataset => ByRow(==(dataset)), :method => ByRow(==(:mle)))
    df_bidra = subset(all_df, :dataset => ByRow(==(dataset)), :method => ByRow(!=(:mle)))

    method_map = :method => renamer([:mle => "Levenberg-\nMarquardt", :bidra_quantile => "BiDRA posteriors\nper quantile"]) => "Methods"
    sub_map = :sub => renamer([:all => "All pairs", :both_c => "Complete pairs\nSD ≥ 20", :both_i => "Incomplete pairs\nSD < 20"])
    row_map = show_row ?
        (:metric => renamer([:LDR => "LDR", :HDR => "HDR", :ic50 => "IC50", :slope => "Slope"])) :
        (:metric => renamer([:LDR => "", :HDR => "", :ic50 => "", :slope => ""]))

    base_mapping = mapping(
        method_map,
        :r_swap,
        row = row_map,
        col = sub_map,
    )

    spec_bar_bidra = data(df_bidra) * base_mapping *
        mapping(dodge = :q => sorter(q_df.sym), color = :q) *
        visual(BarPlot; dodge_gap = 0.001)

    spec_bar_mle = data(df_mle) * base_mapping * visual(BarPlot; gap=1.5, color=HSL(0, 0, 0.3))
    df_median = subset(all_df,
        :dataset => ByRow(==(dataset)),
        :q => ByRow(==(Symbol("0.5"))),
        :method => ByRow(==(:bidra_quantile))
    )
    # spec_hl = data((pos = [0.5, 0.75],)) * mapping(:pos) * visual(HLines; linestyle=:dash)

    spec_hl = data(df_median) * mapping(:r_swap, row=row_map, col=sub_map) * visual(HLines; linestyle=:dash)
    spec_hl += data((pos = [0.],)) * mapping(:pos) * visual(HLines; linestyle=:solid)

    return spec_hl + spec_bar_mle + spec_bar_bidra
end

edge_weight(x, γ=1; mid=0.0, edge=1.0) = mid + (edge - mid) * (2abs(x - 0.5))^γ

function my_col(x)
    s = edge_weight(x, 0.5; mid=0., edge=1.0)
    h = x < 0.5 ? 240. : 0.
    l = edge_weight(x, 1; mid=0.0, edge=0.9)
    return HSL(h, s, l)
end    

set_theme!(Theme(
    fonts = (
        regular = "Noto Sans",
        bold = "Noto Sans Bold"
)))

function prep_figure(dts, dict_all_df, cor_fn, q_df)
    fig = Figure(size=(2200, 1800))

    specs = [dataset_spec(dict_all_df[cor_fn], dts[i, :sym], show_row = (i==nrow(dts)), q_df) for i ∈ 1:nrow(dts)]

    for (i, spec) ∈ enumerate(specs)
        axis = (
            width = 200, height = 200,
            xticklabelrotation = pi / 4,
            xticklabelsize = 18.,
            xgridvisible = false,
            xlabel = "",
            yminorticks = -0.1:0.1:1.0,
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
                yminorticksvisible = true,
                ylabel = "Orientation-averaged $cor_fn correlation coefficient",
                ylabelsize = 24
                )
            else
                axis = (axis...,
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

    resize_to_layout!(fig)
    return fig
end

function AlgebraOfGraphics.aesthetic_mapping(::Type{<:Hexbin},
    ::AlgebraOfGraphics.Normal,
    ::AlgebraOfGraphics.Normal)

        return AlgebraOfGraphics.dictionary([
        1 => AlgebraOfGraphics.AesX,
        2 => AlgebraOfGraphics.AesY,
    ])
end