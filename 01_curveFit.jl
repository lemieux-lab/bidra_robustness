using DataFrames
using CSV
using Optim, GLM, LsqFit

include("utils.jl")

###### Global Var ####
juliaMLe_fn = "/home/golem/scratch/labellec/_RESULTS/MLE_ESTIMATES/all_julia_curveFit.csv"
@. model(x, p) = p[2] + (p[1] - p[2]) / (1 + exp(p[4] * (x - p[3])))
p₀ = [100.,0.,0.,1]

###### Data ##########
data_prefix = "/home/golem/scratch/labellec/_DATA/"
data_df = getRawData(["gCSI", "gray", "ctrpv2"], data_prefix, true)
exp_id = unique(data_df.exp_id)

for e in exp_id
    println(e)
    e_tmp = filter(:exp_id => x -> x == e, data_df)
    e_dataset = unique(e_tmp.dataset)[1]

    xdata = e_tmp.Concentration
    ydata = e_tmp.Viability
    fit = curve_fit(model, xdata, ydata, p₀)

    bestFitParam = fit.param
    rmsd = sqrt(sum(fit.resid .^ 2) / length(fit.resid))
    convergence = fit.converged
    aac = computeAAC(xdata, ydata, bestFitParam)

    tmp_df = DataFrame(LDR=bestFitParam[1], HDR=bestFitParam[2], ic50=bestFitParam[3],
                       slope=bestFitParam[4], aac=aac, exp_id=e, dataset=e_dataset, rmsd=rmsd, convergence=convergence)

    CSV.write(juliaMLe_fn, tmp_df, delim=",", append=true)

end