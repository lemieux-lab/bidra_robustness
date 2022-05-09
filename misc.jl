using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig

include("/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/utils.jl")

###
### Measuring deviation of viabilities response
f_llogistic = llogistic([100, 0, 2., 3]) 

x1 = [0:0.5:4;]
x2 = [-2:0.5:2.5;]
x3 = [-2.5:0.5:2.;]
x4 = [-3:0.5:1.5;]

y1 = f_llogistic.(x1)
y2 = f_llogistic.(x2)
y3 = f_llogistic.(x3)
y4 = f_llogistic.(x4)

Gadfly.set_default_plot_size(3inch, 2inch)
p1 = Gadfly.plot(x=x1, y=y1, Geom.point(), Coord.Cartesian(ymin=0, ymax=100))
p2 = Gadfly.plot(x=x2, y=y2, Geom.point(), Coord.Cartesian(ymin=0, ymax=100))
p3 = Gadfly.plot(x=x3, y=y3, Geom.point(), Coord.Cartesian(ymin=0, ymax=100))
p4 = Gadfly.plot(x=x4, y=y4, Geom.point(), Coord.Cartesian(ymin=0, ymax=100))

std.([y1, y2, y3, y4])
mad.([y1, y2, y3, y4])
maximum.([y1, y2, y3, y4]) .- minimum.([y1, y2, y3, y4])

### Measuring deviation of posterior
N = 4000

x1 = rand(Normal(0,10), N)
x2 = rand(Normal(100,10), N)
x3 = rand(MixtureModel([SkewNormal(0, 10, 1), 
                        Uniform(0, 100), 
                        SkewNormal(100, 20, -5)], [0.4, 0.5, 0.1]), N)
x4 = rand(LogNormal(0.5, 1), N)

Gadfly.set_default_plot_size(3inch, 2inch)
p1 = Gadfly.plot(x=x1, Geom.histogram())
p2 = Gadfly.plot(x=x2, Geom.histogram())
p3 = Gadfly.plot(x=x3, Geom.histogram())
p4 = Gadfly.plot(x=x4, Geom.histogram())

std.([x1, x2, x3, x4])
mad.([x1, x2, x3, x4])
percentile.([x1, x2, x3, x4], 97.5) .- percentile.([x1, x2, x3, x4], 2.5)

#### Testing AAC calculation
expIdSubset_list = ["NCI-H1648_AZ-628_8h"]

data_df = getRawData(["gCSI"], "/home/golem/scratch/labellec/_DATA/", true)
data_df = filter(:exp_id => x -> x ∈ expIdSubset_list, data_df)

ml_df = getMLestimates(["gCSI"], false, missing)
ml_df = filter(:exp_id => i -> i ∈ expIdSubset_list, ml_df)
computeAAC(data_df.Concentration, data_df.Viability, ml_df[1, [:LDR, :HDR, :ic50, :slope]])

posterior_df = getBIDRAposterior("gCSI", expIdSubset_list)
aacPosterior = []
for i in 1:nrow(posterior_df)
    tmp = computeAAC(data_df.Concentration, data_df.Viability, posterior_df[i, [:LDR, :HDR, :ic50, :slope]])
    push!(aacPosterior, tmp)
end
