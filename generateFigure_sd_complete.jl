using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig
using CairoMakie

include("utils.jl")

figure_prefix = "_generated_figures/supp_fig/sd_complete/";

dt = "ctrpv2"

@time expId_list = getExpId_h5(dt);
@time si = StrIndex(expId_list);
@time data_df = getRawData_h5(dt, false, si);

sd_df = combine(groupby(data_df, :exp_id), :Viability => std => :sd) 
expId_complete = filter(:sd => x -> x ≥ 20, sd_df)[:, :exp_id]

fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1, 1], xlabel="Response SD (%)")
CairoMakie.hist!(ax, sd_df.sd, bins=25, color=:black)
CairoMakie.vlines!(ax, [20], color=:blue)
CairoMakie.text!(ax, 0.6, 0.8, text="N = $(nrow(sd_df))\nn(complete) = $(length(expId_complete))", color=:blue, space=:relative,)
fn = "$(figure_prefix)$(dt)_sd_distribution.pdf"
CairoMakie.save(fn, fig)

response_complete = filter(:exp_id => x -> x ∈ expId_complete, data_df)
response_incomplete = filter(:exp_id => x -> x ∉ expId_complete, data_df)
sort!(response_complete, :Concentration)
sort!(response_incomplete, :Concentration)

function plotViabilityByDoseRank(response_df::DataFrame)
    nDose_by_expId = combine(groupby(response_df, :exp_id), nrow=> :nDose)
    nDose = unique(nDose_by_expId.nDose)

    fig = CairoMakie.Figure()

    for j in 1:length(nDose)
        n = nDose[j]
        ax = CairoMakie.Axis(fig[j, 1])

        tmp_expId = filter(:nDose => x -> x == n, nDose_by_expId)[:, :exp_id]
        tmp_n = filter(:exp_id => x -> x ∈ tmp_expId, response_df)

        CairoMakie.boxplot!(ax, tmp_n.conc_id, tmp_n.Viability)
        CairoMakie.ylims!(ax, [0,200])
        #for i in 1:length(unique(tmp_n.conc_id))
        #   tmp = filter(:conc_id => x -> x == i, tmp_n)
        #    CairoMakie.hist!(ax, tmp.Viability, scale_to=-0.6, offset=i, direction=:x)
        #end
    end

    return fig
end

### Complete response viability 
tmp = groupby(response_complete, :exp_id)
for gp in tmp
    gp[!, :conc_id] = StatsBase.denserank(gp.Concentration)
end
response_complete = DataFrame(tmp)
fig = plotViabilityByDoseRank(response_complete);

fn = "$(figure_prefix)$(dt)_complete_exp_response_by_conc_rank.pdf"
CairoMakie.save(fn, fig)

### Incomplete response viability 
tmp = groupby(response_incomplete, :exp_id)
for gp in tmp
    gp[!, :conc_id] = StatsBase.denserank(gp.Concentration)
end
response_incomplete = DataFrame(tmp)
fig = plotViabilityByDoseRank(response_incomplete);

fn = "$(figure_prefix)$(dt)_incomplete_exp_response_by_conc_rank.pdf"
CairoMakie.save(fn, fig)



#### Completeness correlation
pairings_df = getPairings_h5(dt, si);
pairings_df.pairID = collect(1:nrow(pairings_df));

rep1 = innerjoin(pairings_df[:, [:rep_1, :pairID]], sd_df, on=:rep_1=>:exp_id)
rep2 = innerjoin(pairings_df[:, [:rep_2, :pairID]], sd_df, on=:rep_2=>:exp_id)
pair_completeness = innerjoin(rep1, rep2, on=:pairID, makeunique=true)

fig = CairoMakie.Figure();
ax = CairoMakie.Axis(fig[1, 1]);
CairoMakie.hspan!(ax, 20, 60, color=(:pink, 0.2))
CairoMakie.vspan!(ax, 20, 60, color=(:pink, 0.2))
CairoMakie.hexbin!(ax, pair_completeness.sd, pair_completeness.sd_1, bins=75, colormap=["gray95", "gray35", "black"])
CairoMakie.ablines!(ax, 0, 1)

CairoMakie.ylims!(ax, [0,60])
CairoMakie.xlims!(ax, [0,60])

fn = "$(figure_prefix)$(dt)_sd_comparison_biological_replicates.pdf"
CairoMakie.save(fn, fig)