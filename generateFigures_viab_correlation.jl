using DataFrames
using Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig
using ProgressBars
using JuBox

include("utils.jl")
fn = "_generated_figures/supp_fig/viab_corr/"
scaleColor = Scale.lab_gradient("gray95", "black")

### Define dataset to analyze
#dt = ARGS[1]
dt = "ctrpv2"

### Get data
println("1. Get data")
expId_list = getExpId_h5(dt);
si = StrIndex(expId_list);
pairings_df = getPairings_h5(dt, si);
pairings_df.pairID = collect(1:nrow(pairings_df));
viability_df = getRawData_h5(dt, false, si);

println("2. Pair experiments responses")
rep1 = innerjoin(pairings_df[:, [:rep_1, :pairID]], viability_df, on=:rep_1=>:exp_id)
rep2 = innerjoin(pairings_df[:, [:rep_2, :pairID]], viability_df, on=:rep_2=>:exp_id)
pair_shared_dose = innerjoin(rep1, rep2, on=[:pairID, :Concentration], makeunique=true)

function rmse(res)
    return sqrt(sum(res .^ 2) / length(res))
end

if dt ∈ ["gCSI", "ctrpv2"]
    println("1. Plot RMSΔ by # of concentration")
    local conc_count = combine(groupby(pair_shared_dose, :pairID), :Concentration => length∘unique => :nConcentration)
    diff_conc = combine(groupby(pair_shared_dose, [:pairID, :Concentration]), [:Viability, :Viability_1] => ((a, b) -> a-b) => :Δ)
    local rmsΔ = combine(groupby(diff_conc, :pairID), :Δ => rmse => :rmsΔ)

    rmsΔ = innerjoin(conc_count, rmsΔ, on=:pairID)
    rmsΔ_complete = filter(:nConcentration => x -> x == maximum(rmsΔ.nConcentration), rmsΔ)

    Gadfly.set_default_plot_size(4inch, 3inch)
    local p_rmsΔ = Gadfly.plot(rmsΔ, x=:nConcentration, y=:rmsΔ, yintercept=[median(rmsΔ.rmsΔ), median(rmsΔ_complete.rmsΔ)],
                    Geom.boxplot(), Geom.hline(),
                    Scale.x_discrete)
    println(countmap(rmsΔ.nConcentration))
    draw(PDF(fn*"viability_rmse_$dt.pdf", 4inch, 3inch), p_rmsΔ)

    println("2. Plot viability rep. correlation")
    Gadfly.set_default_plot_size(3inch, 2inch)
    p_repViab = Gadfly.plot(pair_shared_dose, x=:Viability, y=:Viability_1,
                 Geom.abline(style=:dash), Geom.hexbin(xbincount=80, ybincount=80),
                 #Scale.color_continuous(colormap=scaleColor, minvalue=1),
                 Scale.color_log10(colormap=scaleColor, minvalue=1),
                 Coord.cartesian(xmin=-50, ymin=-50, xmax=200, ymax=200),
                 Theme(panel_stroke="black"),
                 Guide.title(dt));
    draw(PDF(fn*"viability_corr_$dt.pdf", 4inch, 3inch), p_repViab)
else 
    R = 10000
    rep_rmsΔ_df = DataFrame()

    local conc_count = combine(groupby(pair_shared_dose, :pairID), :Concentration => length∘unique => :nConcentration)

    println("1. Run bootstrap for R=$R")
    for r in ProgressBar(1:R)
        tmp = mapreduce(g -> DataFrame(g[rand(1:nrow(g)), :]), vcat, groupby(pair_shared_dose, [:pairID, :Concentration]))
        diff_conc = combine(groupby(tmp, [:pairID, :Concentration]), [:Viability, :Viability_1] => ((a, b) -> a-b) => :Δ)
        local rmsΔ = combine(groupby(diff_conc, :pairID), :Δ => rmse => :rmsΔ)
        rmsΔ = innerjoin(conc_count, rmsΔ, on=:pairID)
        append!(rep_rmsΔ_df, rmsΔ)

        if r == 1 
            println("--> Plot viability rep. correlation for random sampling ")
            Gadfly.set_default_plot_size(3inch, 2inch)
            local p_repViab = Gadfly.plot(tmp, x=:Viability, y=:Viability_1,
                        Geom.abline(style=:dash), Geom.hexbin(xbincount=80, ybincount=80),
                        #Scale.color_continuous(colormap=scaleColor, minvalue=1),
                        Scale.color_log10(colormap=scaleColor, minvalue=1),
                        Coord.cartesian(xmin=-50, ymin=-50, xmax=200, ymax=200),
                        Theme(panel_stroke="black"),
                        Guide.title(dt));
            draw(PDF(fn*"viability_corr_$dt.pdf", 4inch, 3inch), p_repViab)
        end
    end

    println("3. Combine all sampling")
    rep_mean_df = combine(groupby(rep_rmsΔ_df, :pairID), :nConcentration => unique => :nConcentration, :rmsΔ => mean => :rmsΔ_mean, :rmsΔ => (x -> percentile(x, 2.5,)) => :low, :rmsΔ => (x -> percentile(x, 97.5,)) => :up, :rmsΔ => std => :std)
    sort!(rep_mean_df, :rmsΔ_mean)

    Gadfly.set_default_plot_size(6inch, 3inch)
    #p = Gadfly.plot(rep_mean_df, x=:pairID, y=:rmsΔ_mean, ymin=:low, ymax=:up, color=:nConcentration,
    #                Geom.errorbar(), Geom.bar(),
    #                Scale.x_discrete);

    println("4. Plot RMSΔ by # of concentration")
    local p_rmsΔ = Gadfly.plot(rep_mean_df, x=:rmsΔ_mean, y=:std, xgroup=:nConcentration, Geom.subplot_grid(Geom.point()));

    println(countmap(rep_mean_df.nConcentration))
    draw(PDF(fn*"viability_rmse_$dt.pdf", 6inch, 3inch), p_rmsΔ)
end

