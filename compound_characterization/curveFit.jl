using DataFrames
using CSV
using Optim, GLM, LsqFit

include("../utils.jl")

###### Global Var ####
juliaMLe_fn = "../data/all_julia_curveFit.csv"
@. model(x, p) = p[2] + ((p[1] - p[2]) / (1 + 10^(p[4] * (x - p[3]))))
p₀ = [100.,0.,0.,1]

###### Data ##########
for dt in ["gCSI", "gray", "ctrpv2"]
    data_df = getRawData_h5(dt, false)
    exp_id = getExpId_h5(dt)

    for e in exp_id
        println(e)
        e_tmp = filter(:exp_id => x -> x == e, data_df)

        xdata = e_tmp.Concentration
        ydata = e_tmp.Viability
        fit = curve_fit(model, xdata, ydata, p₀)

        bestFitParam = fit.param
        rmse = sqrt(sum(fit.resid .^ 2) / length(fit.resid))
        convergence = fit.converged
        aac = computeAAC(xdata, ydata, bestFitParam)

        tmp_df = DataFrame(LDR=bestFitParam[1], HDR=bestFitParam[2], ic50=bestFitParam[3],
                        slope=bestFitParam[4], aac=aac, exp_id=e, dataset=dt, rmse=rmse, convergence=convergence)

        CSV.write(juliaMLe_fn, tmp_df, delim=",", append=true)
    end
end