using JuBox
include("../utils.jl")
include("aac_plot.jl")

figure_prefix = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_figures/supp_fig/aac_analysis/"

### Estimates and infer AAC for common range for pairs of biological replicates
dt = "ctrpv2"
expId_list = getExpId_h5(dt)
si = StrIndex(expId_list)
data_df = getRawData_h5(dt, false, si)
ml_df = getMLestimates(dt, si)
posterior_df = getPosterior_h5(dt, false, si)

pairings_df = getPairings_h5(dt, si)
m = nrow(pairings_df)

ind_range = combine(groupby(data_df, :exp_id), :Concentration => minimum => :min,
                                               :Concentration => maximum => :max)
#sort!(ind_range, [:min, :max])

## Get shared range by pairs
paired_min = Array{Float64, 1}(undef, m)
paired_max = Array{Float64, 1}(undef, m)

for i in 1:m 
    rep1, rep2 = Array(pairings_df[i, :])
    paired_min[i] = maximum([ind_range[ind_range.exp_id .== rep1, :min][1], ind_range[ind_range.exp_id .== rep2, :min][1]])
    paired_max[i] = minimum([ind_range[ind_range.exp_id .== rep1, :max][1], ind_range[ind_range.exp_id .== rep2, :max][1]])
end

pairings_df[!, :min] = paired_min
pairings_df[!, :max] = paired_max

## plot pairing range
f = plot_concentrationRangeOverview(pairings_df, ind_range);
fn = "$(figure_prefix)min_max_concentration_range_$dt.pdf"
CairoMakie.save(fn, f)

## Calculate AAC
tmp1 = innerjoin(ml_df, pairings_df[:, [:rep_1, :min, :max]], on=:exp_id=>:rep_1)
tmp2 = innerjoin(ml_df, pairings_df[:, [:rep_2, :min, :max]], on=:exp_id=>:rep_2)
ml_df = vcat(tmp1, tmp2)
ml_df[:, :aacRC] = map(r -> computeAAC(r[10:11], r[1:4]), eachrow(ml_df))
ml_df[:, :Δaac] = abs.(ml_df.aac .- ml_df.aacRC)

tmp1 = innerjoin(posterior_df, pairings_df[:, [:rep_1, :min, :max]], on=:exp_id=>:rep_1)
tmp2 = innerjoin(posterior_df, pairings_df[:, [:rep_2, :min, :max]], on=:exp_id=>:rep_2)
posterior_df = vcat(tmp1, tmp2)
posterior_df[:, :aacRC] = map(r -> computeAAC(r[9:10], r[1:4]), eachrow(posterior_df))

median_df = combine(groupby(posterior_df, :exp_id), :aacRC => median => :aacMedian)

## Look for Inf and NAN
ml_not = filter(:aacRC => x -> !isfinite(x), ml_df)
bidra_not = filter(:aacRC => x -> !isfinite(x), posterior_df)
median_not = filter(:aacMedian => x -> !isfinite(x), median_df)

f = plot_experimentWithNoAAC(ml_not, data_df);
fn = "$(figure_prefix)lev_marq_aac_nan_$dt.pdf"
CairoMakie.save(fn, f)

#### Compare AAC
sd_viab = combine(groupby(data_df, :exp_id), :Viability => std => :viab_sd)

## Levenberg-Marquardt
paired_ml = innerjoin(pairings_df, ml_df[:, [:exp_id, :aacRC]], on=:rep_1 => :exp_id, renamecols = ""=>"_1")
paired_ml = innerjoin(paired_ml, ml_df[:, [:exp_id, :aacRC]], on=:rep_2 => :exp_id, renamecols = ""=>"_2")
paired_ml = filter([:aacRC_1, :aacRC_2] => (a, b) -> isfinite(a) && isfinite(b), paired_ml)
paired_ml = innerjoin(paired_ml, sd_viab, on=:rep_1=>:exp_id, renamecols=""=>"_1")
paired_ml = innerjoin(paired_ml, sd_viab, on=:rep_2=>:exp_id, renamecols=""=>"_2")

## Median
paired_median = innerjoin(pairings_df, median_df[:, [:exp_id, :aacMedian]], on=:rep_1 => :exp_id, renamecols = ""=>"_1")
paired_median = innerjoin(paired_median, median_df[:, [:exp_id, :aacMedian]], on=:rep_2 => :exp_id, renamecols = ""=>"_2")
paired_median = filter([:aacMedian_1, :aacMedian_2] => (a, b) -> isfinite(a) && isfinite(b), paired_median)
paired_median = innerjoin(paired_median, sd_viab, on=:rep_1=>:exp_id, renamecols=""=>"_1")
paired_median = innerjoin(paired_median, sd_viab, on=:rep_2=>:exp_id, renamecols=""=>"_2")

## posterior
posterior1 = mapreduce(e -> filter(:exp_id => x -> x == e, posterior_df), vcat, pairings_df.rep_1)
posterior1 = innerjoin(posterior1, sd_viab, on=:exp_id=>:exp_id, renamecols="_1"=>"_1")

posterior2 = mapreduce(e -> filter(:exp_id => x -> x == e, posterior_df), vcat, pairings_df.rep_2)
posterior2 = innerjoin(posterior2, sd_viab, on=:exp_id=>:exp_id, renamecols="_2"=>"_2")

paired_posterior = hcat(posterior1[:, [:exp_id, :aacRC_1, :viab_sd_1]], posterior2[:, [:exp_id, :aacRC_2, :viab_sd_2]], makeunique=true)
rename!(paired_posterior, :exp_id => :rep_1)
rename!(paired_posterior, :exp_id_1 => :rep_2)

posteriorQQ = combine(groupby(paired_posterior, [:rep_1, :rep_2]), :aacRC_1 => sort => :sorted_1, 
                                                                   :aacRC_2 => sort => :sorted_2,
                                                                   :viab_sd_1 => :viab_sd_1,
                                                                   :viab_sd_2 => :viab_sd_2)

f = plot_AACcomparison(paired_ml, paired_median, posteriorQQ);
fn = "$(figure_prefix)aac_corr_all_$dt.pdf"
CairoMakie.save(fn, f)

#### Explore negative AAC for posterior
## Only considering experiments that are part of a pair of replicates
negAAC_posterior = filter(:aacRC => x -> x < 0, posterior_df)
negAAC_posterior = innerjoin(negAAC_posterior, sd_viab, on=:exp_id)
sample_expId = combine(groupby(negAAC_posterior, :exp_id), :aacRC => length => :nSample,
                                                           :viab_sd => unique => :viab_sd)
sample_expId[!, :group] = map(sd -> sd < 20 ? 2 : 1, sample_expId.viab_sd)
sample_expId[!, :propIte] = sample_expId.nSample ./ 4000

f = plot_negativeAACposterior(negAAC_posterior, sample_expId);
fn = "$(figure_prefix)neg_aac_$dt.pdf"
CairoMakie.save(fn, f)

