#push!(LOAD_PATH, "Utils/")
using Statistics, StatsBase
using Distributions
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("utils.jl")

###### Global Var ####
figure_prefix = "_generated_figures/model/"
datasets= ["gCSI"]#["gray", "gCSI", "ctrpv2"]
expId_list = mapreduce(dt -> getExpId_h5(dt), vcat, datasets);
nId_list = map(dt -> length(getExpId_h5(dt)), datasets);
si = StrIndex(expId_list)

###### Data ##########
data_df = mapreduce(dt -> getRawData_h5(dt, false, si), vcat, datasets);
ml_df = mapreduce(dt -> getMLestimates(dt, si), vcat, datasets);

### example experiments from gCSI
expId_subset_list = ["NCI-H1648_AZ-628_8h", "Calu-1_PF-4708671_6b", "RERF-LC-MS_Gemcitabine_4b", "HCC78_Lapatinib_11a"];
data_subset = filter(:exp_id => x -> si.id2str[x] ∈ expId_subset_list, data_df);
ml_subset = filter(:exp_id => x -> si.id2str[x] ∈ expId_subset_list, ml_df);
posterior_subset = getPosterior_h5("gCSI", false, si, expId_subset_list);

### Raw data overlook
maxConcentration_df = combine(groupby(data_df, :exp_id), 
        :Concentration => maximum => :Concentration, 
        [:Concentration, :Viability] => ((a,b) -> b[argmax(a)])  => :Viability,
        [:Concentration, :dataset] => ((a,b) -> b[argmax(a)])  => :dataset);

minConcentration_df = combine(groupby(data_df, :exp_id), 
        :Concentration => minimum => :Concentration, 
        [:Concentration, :Viability] => ((a,b) -> b[argmin(a)])  => :Viability,
        [:Concentration, :dataset] => ((a,b) -> b[argmin(a)])  => :dataset)

## Draw Viability of max [] distributions per dataset
Gadfly.set_default_plot_size(4inch, 3inch)
x = percentile(maxConcentration_df.Viability, [2.5, 97.5])
ymin = [0.,0.]
ymax = ymin .+ 0.20
p1 = Gadfly.plot(maxConcentration_df, x=:Viability, color=:dataset, 
                Geom.histogram(position=:stack, bincount=200, density=true),
                layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon()), 
                Coord.cartesian(xmin=-20, xmax=150),
                Theme(panel_stroke="black"))
draw(PDF(figure_prefix*"viability_max_conc_LARGE.pdf", 4inch, 3inch), p1)

## Draw Viability of min [] distributions per dataset
Gadfly.set_default_plot_size(4inch, 3inch)
x = percentile(minConcentration_df.Viability, [2.5, 97.5])
ymin = [0.,0.]
ymax = ymin .+ 0.20
p2 = Gadfly.plot(minConcentration_df, x=:Viability, color=:dataset, 
                Geom.histogram(position=:stack, bincount=200, density=true),
                layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon()),
                Coord.cartesian(xmin=0, xmax=150),
                Theme(panel_stroke="black"))
draw(PDF(figure_prefix*"viability_min_conc_LARGE.pdf", 4inch, 3inch), p2)

## Draw [] distributions per dataset
Gadfly.set_default_plot_size(4inch, 3inch)
x = percentile(data_df.Concentration, [2.5, 97.5])
ymin = [0.,0.]
ymax = ymin .+ 2.
p3 = Gadfly.plot(data_df, x=:Concentration, color=:dataset,
                Geom.histogram(position=:stack, density=true),
                layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon()),
                Coord.cartesian(xmin=-7),
                Theme(panel_stroke="black"))
draw(PDF(figure_prefix*"concentration_LARGE.pdf", 4inch, 3inch), p3)


### Prior density plots
N = 400000
## IC50
ic50_prior = rand(Normal(0,10), N)
x = percentile(data_df.Concentration, [2.5, 97.5])
ymin = [0.,0.]
ymax = ymin .+ 0.04
p4 = Gadfly.plot(layer(x=ic50_prior, Geom.density(bandwidth=1)),
                 layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon(),
                 Theme(panel_stroke="black")))
draw(PDF(figure_prefix*"ic50_prior_LARGE.pdf", 4inch, 3inch), p4)

## LDR
ldr_prior = rand(Normal(100,10), N)
x = percentile(minConcentration_df.Viability, [2.5, 97.5])
ymin = [0.,0.]
ymax = ymin .+ 0.04
p5 = Gadfly.plot(layer(x=ldr_prior, Geom.density(bandwidth=1)),
                 layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon(),
                 ))
draw(PDF(figure_prefix*"ldr_prior.pdf", 4inch, 3inch), p5)

## HDR
λₖ = [0.4, 0.5, 0.1]
hdr_mm = MixtureModel([SkewNormal(0, 10, 1), Uniform(0, 100), SkewNormal(100, 20, -5)], λₖ)
hdr_prior = rand(hdr_mm, N)
x = percentile(maxConcentration_df.Viability, [2.5, 97.5])
ymin = [0.,0.]
ymax = ymin .+ 0.04
p6 = Gadfly.plot(layer(x=hdr_prior, Geom.density(bandwidth=1)),
                 layer(x=x, ymin=ymin, ymax=ymax, Geom.ribbon() ))
draw(PDF(figure_prefix*"hdr_prior.pdf", 4inch, 3inch), p6)

## slope
slope_prior = rand(LogNormal(0.5, 1), N)
p7 = Gadfly.plot(layer(x=slope_prior, Geom.density(bandwidth=1)),
                Coord.cartesian(xmin=-5, xmax=40))
draw(PDF(figure_prefix*"slope_prior.pdf", 4inch, 3inch), p7)

## σ
σ_prior = rand(LogNormal(1, 1), N)
p8 = Gadfly.plot(layer(x=σ_prior, Geom.density(bandwidth=1)),
                Coord.cartesian(xmin=-5, xmax=50))
draw(PDF(figure_prefix*"sigma_prior.pdf", 4inch, 3inch), p8)


### Examples of analysis
Gadfly.set_default_plot_size(6inch, 4inch)
p9 = Gadfly.plot(layer(data_subset, x=:Concentration, y=:Viability, color=si.id2str[data_subset.exp_id], Geom.point()), Scale.color_discrete())

for e in expId_subset_list
        println(e)
        postCurves = getPosteriorCurves(filter(:exp_id => x -> si.id2str[x] == e, posterior_subset), -4, 1)
        xDose = parse.(Float64, names(postCurves))

        med_curve = median.(eachcol(postCurves))
        upper_curve = percentile.(eachcol(postCurves), 97.5)
        lower_curve = percentile.(eachcol(postCurves), 2.5)
        push!(p9, layer(x=xDose, y=med_curve, color=[e], linestyle=[:dot], Geom.line()))
        
        mleFit = llogistic(Array(filter(:exp_id => x -> si.id2str[x] == e, ml_subset)[1, [:LDR, :HDR, :ic50, :slope]]))
        yViability = mleFit.(xDose)
        push!(p9, layer(x=xDose, y=yViability, color=[e], Geom.line()))
        push!(p9, layer(x=xDose, ymin=lower_curve, ymax=upper_curve, color=[e], alpha=[0.3], Geom.ribbon()))
end
display(p9)
draw(PDF(figure_prefix*"curves_examples.pdf", 6inch, 4inch), p9)

### posterior distributions
### Add ΔHDR
Gadfly.set_default_plot_size(4inch, 2inch)
for e in expId_subset_list
        tmp = filter(:exp_id => x -> si.id2str[x] == e, posterior_subset)
        est = filter(:exp_id => x -> si.id2str[x] == e, ml_subset)[1, :HDR]
        ΔHDR = diff(percentile(tmp.HDR, [2.5, 97.5]))[1]
        p10 = Gadfly.plot(layer(x=hdr_prior, Geom.density(bandwidth=1)),
                         layer(tmp, x=:HDR, Geom.histogram(position=:stack, bincount=200, density=true)),
                         layer(xintercept=[est], Geom.vline(color=colorant"black")),
                         layer(xintercept=percentile(tmp.HDR, [2.5, 97.5]), Geom.vline()),
                         Coord.Cartesian(xmin=-10, xmax=110),
                         Guide.title("ΔHDR = $ΔHDR"))
        draw(PDF(figure_prefix*"HDR_posterior_prior_"*e*".pdf", 4inch, 2inch), p10)
        #display(p10)
end

for e in expId_subset_list
        tmp = filter(:exp_id => x -> si.id2str[x] == e, posterior_subset)
        est = filter(:exp_id => x -> si.id2str[x] == e, ml_subset)[1, :LDR]
        p10 = Gadfly.plot(layer(x=ldr_prior, Geom.density(bandwidth=1)),
                         layer(tmp, x=:LDR, Geom.histogram(position=:stack, bincount=200, density=true)),
                         layer(xintercept=[est], Geom.vline(color=colorant"black")),
                         layer(xintercept=percentile(tmp.LDR, [2.5, 97.5]), Geom.vline()),
                         Coord.Cartesian(xmin=90, xmax=120))
        draw(PDF(figure_prefix*"LDR_posterior_prior_"*e*".pdf", 4inch, 2inch), p10)
        #display(p10)
end

### Add exp.dose-range
for e in expId_subset_list
        tmp = filter(:exp_id => x -> si.id2str[x] == e, posterior_subset)
        est = filter(:exp_id => x -> si.id2str[x] == e, ml_subset)[1, :ic50]
        data_subset = filter(:exp_id => x -> x == si.str2id[e], data_df)
        p10 = Gadfly.plot(layer(x=ic50_prior, Geom.density(bandwidth=1)),
                         layer(tmp, x=:ic50, Geom.histogram(position=:stack, bincount=200, density=true)),
                         layer(xintercept=[est], Geom.vline(color=colorant"black")),
                         layer(xintercept=percentile(tmp.ic50, [2.5, 97.5]), Geom.vline()),
                         layer(xintercept=[minimum(data_subset.Concentration), maximum(data_subset.Concentration)], Geom.vline(color=colorant"orange")),
                         Coord.Cartesian(xmin=-5, xmax=5))
        draw(PDF(figure_prefix*"ic50_posterior_prior_"*e*".pdf", 4inch, 2inch), p10)
        display(p10)
end