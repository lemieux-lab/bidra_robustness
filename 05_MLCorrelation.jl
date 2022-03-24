using DataFrames
using Statistics, StatsBase

include("utils.jl")

results_prefix = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/"

gCSIpairings_df = getPairings("gCSI")
gCSImlPaired_df = getMLestimates(["gCSI"], true, gCSIpairings_df)
bidra_params = ["LDR", "HDR", "ic50", "slope", "aac"]

for pr in bidra_params
    rep1 = Symbol(pr,"_rep_1")
    rep2 = Symbol(pr,"_rep_2")

    mlCorr_df = correlationAnalysis(gCSImlPaired_df[:,rep1], gCSImlPaired_df[:,rep2])
    mlCorr_df[:, :N] = [nrow(gCSImlPaired_df)]
    mlCorr_df[:, :dataset] = ["gCSI"]
    mlCorr_df[:, :param] = [pr]

    CSV.write(results_prefix*"mlCorrelations.csv", mlCorr_df, delim=",", append=true)
end

Gadfly.set_default_plot_size(4inch, 3inch)
scaleColor = Scale.lab_gradient("gray95","black")
pr="ic50"
p = Gadfly.plot(gCSImlPaired_df, x=Symbol(pr,"_rep_1"), y=Symbol(pr,"_rep_1"),
            Geom.hexbin(xbincount=75, ybincount=75), 
            Scale.color_continuous(colormap=scaleColor, minvalue=1))
display(p)

p = Gadfly.plot(gCSIml_df, x=:LDR, )