using DataFrames, CSV
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots, CairoMakie
using Cairo, Fontconfig
using ProgressBars

include("utils_iric.jl")

data_df = getViability()
posterior_path = "results/"
figure_path = "figures/"

data_groups = groupby(data_df, [:anonym_id])
global plot_curves_group = DataFrame()
global viab_data_group = DataFrame()

global plot_by_row = 10
global r = 0
global c = 0
n = length(data_groups)

println(r, c)

for expId in ProgressBar(1:n)
    #println("************", expId)
    global c = c+1 ≤ plot_by_row ? c+1 : 1
    global r = c == 1 ? r+1 : r

    tmp = DataFrame(data_groups[expId])
    tmp_chains = run_BiDRA(tmp)
    comp = unique(tmp.anonym_id)[1]
    fn = string(comp, ".csv")
    CSV.write(joinpath(posterior_path, fn), tmp_chains)

    x_concentration = collect(minimum(tmp.log10)-2:0.1:maximum(tmp.log10)+2)
    posterior_curves = Array{Float64}(undef, 0, length(x_concentration))
    for r in 1:nrow(tmp_chains)
        row_curve = llogistic([tmp_chains.LDR[r], tmp_chains.HDR[r], tmp_chains.ic50[r], tmp_chains.slope[r]])
        y_infer = row_curve.(x_concentration)
        posterior_curves = vcat(posterior_curves, y_infer')
    end

    plot_curves = DataFrame(concentration=x_concentration, med_curve=median.(eachcol(posterior_curves)), up_curve=percentile.(eachcol(posterior_curves), [97.5]), lo_curve=percentile.(eachcol(posterior_curves), [2.5]), compound=repeat([comp], length(x_concentration)), row=repeat([r], length(x_concentration)), col=repeat([c], length(x_concentration)))

    global plot_curves_group = vcat(plot_curves_group, plot_curves)
    global viab_data_group = vcat(viab_data_group, tmp)
end

plot_lyt = combine(groupby(plot_curves_group, :compound), :row => unique => :row, :col => unique => :col)
plot_viab = innerjoin(viab_data_group, plot_lyt, on=:anonym_id=>:compound, makeunique=true)

w = plot_by_row
h = maximum(plot_curves_group.row)
sz = 300

f = Figure(backgroundcolor="transparent", resolution=(sz*w, sz*h))
for e in unique(plot_lyt.compound)
    r = plot_lyt[plot_lyt.compound .== e, :row][1]
    c = plot_lyt[plot_lyt.compound .== e, :col][1]
    tmp_viab = filter(:anonym_id => x -> x == e, plot_viab)
    tmp_curves =filter(:compound => x -> x == e, plot_curves_group)

    ax = Axis(f[r, c], title=e)
    CairoMakie.band!(ax, tmp_curves.concentration, tmp_curves.lo_curve, tmp_curves.up_curve, color=(:black, 0.3))
    CairoMakie.scatter!(ax, tmp_viab.log10, tmp_viab.inhib_col, color=:black)
    CairoMakie.lines!(ax, tmp_curves.concentration, tmp_curves.med_curve, color=:black)

       CairoMakie.ylims!(ax, -10, 110)
end
    
fn = string("iric_dr_curves.pdf")
CairoMakie.save("$figure_path$fn", f)