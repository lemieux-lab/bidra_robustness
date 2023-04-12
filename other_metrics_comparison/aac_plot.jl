using CairoMakie

function plot_concentrationRangeOverview(pairings_df::DataFrame, ind_range::DataFrame)
    f = Figure(backgroundcolor="transparent", resolution = (1200, 800));
    axMin = Axis(f[1, 1], title="Individual Min. Concentration", xlabel="Rep. 1", ylabel="Rep. 2")
    axMax = Axis(f[1, 2], title="Individual Max. Concentration", xlabel="Rep. 1", ylabel="Rep. 2")
    axCountMin = Axis(f[3, 1], title="Min. Concentration", xlabel="Concentration")
    axCountMax = Axis(f[3, 2], title="Max. Concentration", xlabel="Concentration")
    axRange = Axis(f[1:3, 3], title="Shared Concentration range", xlabel="Concentration", ylabel="Exp. Pair")

    rep1_val = filter(:exp_id => x -> x ∈ pairings_df.rep_1, ind_range)
    rep2_val = filter(:exp_id => x -> x ∈ pairings_df.rep_2, ind_range)

    hb_min = hexbin!(axMin, rep1_val.min, rep2_val.min, cellsize=0.1, threshold=1, colormap=["gray95", "gray35", "black"])
    ablines!(axMin, [0], [1])
    Colorbar(f[2, 1], hb_min, vertical=false, label="Pairs Count")

    hb_max = hexbin!(axMax, rep1_val.max, rep2_val.max, cellsize=0.1, threshold=1, colormap=["gray95", "gray35", "black"])
    ablines!(axMax, [0], [1])
    Colorbar(f[2, 2], hb_max, vertical=false, label="Pairs Count")

    hist!(axCountMin, pairings_df.min, color=:black)
    hist!(axCountMax, pairings_df.max, color=:black)

    pairings_df[!, :Δ] = pairings_df.max .- pairings_df.min
    sortedPairings = sort(pairings_df, [:Δ, :min])
    rangebars!(axRange, 1:m, sortedPairings.min, sortedPairings.max, direction=:x, color=:black, linewidth=0.1)
    return f
end

function plot_experimentWithNoAAC(ml_not::DataFrame, data_df::DataFrame)
    col = 6
    n = nrow(ml_not)
    f = Figure(backgroundcolor="transparent", resolution = (400*col, 400*ceil(n/col)))
    r = 1
    c = 1

    for e in unique(ml_not.exp_id)
        ax = Axis(f[r, c], title=si.id2str[e])

        tmp = filter(:exp_id => x -> x == e, data_df)
        scatter!(ax, tmp.Concentration, tmp.Viability, color=:black)

        curveX = collect(minimum(tmp.Concentration):0.1:maximum(tmp.Concentration))
        param = Array(ml_not[ml_not.exp_id .== e, [:LDR, :HDR, :ic50, :slope, :aacRC]])
        fx = llogistic(param[1:4])
        curveY = fx.(curveX)
        lines!(ax, curveX, curveY, color=:black)
        text!(ax, median(curveX), 25, text="LDR = $(param[1])\nHDR = $(param[2])\nIC50 = $(param[3])\nslope = $(param[4])\naac = $(param[5])", align=[:center, :baseline])

        rangeX = Array(filter(:exp_id => x -> x == e, ml_not)[:, [:min, :max]])
        vlines!(ax, rangeX, color=:black, linestyle=:dash)

        ylims!(ax, [-10, 120])

        c = c == col ? 1 : c+1
        r = c == 1 ? r+1 : r
    end 

    return f
end

function plot_AACcomparison(paired_ml::DataFrame, paired_median::DataFrame, paired_posterior::DataFrame)
    f = Figure(backgroundcolor="transparent", resolution = (450*3, 450*3))
    posX = -50
    posY = 80
    alignArg = [:left, :top]
    cs = 3
    cl_map = ["gray95", "gray35", "black"]

    ### Correlation of Lev.-Marq.
    ## All pairs
    mlCorr_all = correlationAnalysis(paired_ml[:,:aacRC_1], paired_ml[:,:aacRC_2])
    
    ax1 = Axis(f[1, 1], title="Levenberg_Marquardt\nAll pairs(n=$(nrow(paired_ml)))", xlabel="AAC Rep1", ylabel="AAC Rep2")
    hexbin!(ax1, paired_ml.aacRC_1, paired_ml.aacRC_2, cellsize=cs, threshold=1, colormap=cl_map)
    ablines!(ax1, [0], [1])
    text!(ax1, posX, posY, text="rₛ=$(mlCorr_all[1, :rₛ])\nr=$(mlCorr_all[1, :r])", align=alignArg)

    ## Complete pairs
    completePairs_ml = filter([:viab_sd_1, :viab_sd_2] => (a, b) -> a ≥ 20 && b ≥ 20, paired_ml)
    mlCorr_complete = correlationAnalysis(completePairs_ml[:,:aacRC_1], completePairs_ml[:,:aacRC_2])

    ax2 = Axis(f[2, 1], title="Levenberg_Marquardt\nComplete pairs(n=$(nrow(completePairs_ml)))", xlabel="AAC Rep1", ylabel="AAC Rep2")
    hexbin!(ax2, completePairs_ml.aacRC_1, completePairs_ml.aacRC_2, cellsize=cs, threshold=1, colormap=cl_map)
    ablines!(ax2, [0], [1])
    text!(ax2, posX, posY,text="rₛ=$(mlCorr_complete[1, :rₛ])\nr=$(mlCorr_complete[1, :r])", align=alignArg)

    ## incomplete pairs
    incompletePairs_ml = filter([:viab_sd_1, :viab_sd_2] => (a, b) -> a < 20 && b < 20, paired_ml)
    mlCorr_incomplete = correlationAnalysis(incompletePairs_ml[:,:aacRC_1], incompletePairs_ml[:,:aacRC_2])

    ax3 = Axis(f[3, 1], title="Levenberg_Marquardt\nIncomplete pairs(n=$(nrow(incompletePairs_ml)))", xlabel="AAC Rep1", ylabel="AAC Rep2")
    hexbin!(ax3, incompletePairs_ml.aacRC_1, incompletePairs_ml.aacRC_2, cellsize=cs, threshold=1, colormap=cl_map)
    ablines!(ax3, [0], [1])
    text!(ax3, posX, posY, text="rₛ=$(mlCorr_incomplete[1, :rₛ])\nr=$(mlCorr_incomplete[1, :r])", align=alignArg)

    ### Correlation of posterior Median
    ## All pairs
    medianCorr_all = correlationAnalysis(paired_median[:,:aacMedian_1], paired_median[:,:aacMedian_2])
    
    ax4 = Axis(f[1, 2], title="Posterior medians\nAll pairs(n=$(nrow(paired_median)))", xlabel="AAC Rep1", ylabel="AAC Rep2")
    hexbin!(ax4, paired_median.aacMedian_1, paired_median.aacMedian_2, cellsize=cs, threshold=1, colormap=cl_map)
    ablines!(ax4, [0], [1])
    text!(ax4, posX, posY, text="rₛ=$(medianCorr_all[1, :rₛ])\nr=$(medianCorr_all[1, :r])", align=alignArg)

    ## Complete pairs
    completePairs_median = filter([:viab_sd_1, :viab_sd_2] => (a, b) -> a ≥ 20 && b ≥ 20, paired_median)
    medianCorr_complete = correlationAnalysis(completePairs_median[:,:aacMedian_1], completePairs_median[:,:aacMedian_2])
    
    ax5 = Axis(f[2, 2], title="Posterior medians\nComplete pairs(n=$(nrow(completePairs_median)))", xlabel="AAC Rep1", ylabel="AAC Rep2")
    hexbin!(ax5, completePairs_median.aacMedian_1, completePairs_median.aacMedian_2, cellsize=cs, threshold=1, colormap=cl_map)
    ablines!(ax5, [0], [1])
    text!(ax5, posX, posY, text="rₛ=$(medianCorr_complete[1, :rₛ])\nr=$(medianCorr_complete[1, :r])", align=alignArg)

    ## Incomplete pairs
    incompletePairs_median = filter([:viab_sd_1, :viab_sd_2] => (a, b) -> a < 20 && b < 20, paired_median)
    medianCorr_incomplete = correlationAnalysis(incompletePairs_median[:,:aacMedian_1], incompletePairs_median[:,:aacMedian_2])
    
    ax6 = Axis(f[3, 2], title="Posterior medians\nIncomplete pairs(n=$(nrow(incompletePairs_median)))", xlabel="AAC Rep1", ylabel="AAC Rep2")
    hexbin!(ax6, incompletePairs_median.aacMedian_1, incompletePairs_median.aacMedian_2, cellsize=cs, threshold=1, colormap=cl_map)
    ablines!(ax6, [0], [1])
    text!(ax6, posX, posY, text="rₛ=$(medianCorr_incomplete[1, :rₛ])\nr=$(medianCorr_incomplete[1, :r])", align=alignArg)
    
    ### Correlation of posterior Median
    ## All pairs
    qqCorr_all = correlationAnalysis(paired_posterior[:,:sorted_1], paired_posterior[:,:sorted_2])

    ax7 = Axis(f[1, 3], title="Posterior QQ\nAll pairs(n=$(nrow(paired_median)))", xlabel="AAC Rep1", ylabel="AAC Rep2")
    hexbin!(ax7, paired_posterior.sorted_1, paired_posterior.sorted_2, cellsize=cs, threshold=1, colormap=cl_map)
    ablines!(ax7, [0], [1])
    text!(ax7, posX, posY, text="rₛ=$(qqCorr_all[1, :rₛ])\nr=$(qqCorr_all[1, :r])", align=alignArg)

    ## Complete pairs
    completePairs_qq = filter([:viab_sd_1, :viab_sd_2] => (a, b) -> a ≥ 20 && b ≥ 20, paired_posterior)
    qqCorr_complete = correlationAnalysis(completePairs_qq[:,:sorted_1], completePairs_qq[:,:sorted_2])

    ax8 = Axis(f[2, 3], title="Posterior QQ\nComplete pairs(n=$(nrow(completePairs_median)))", xlabel="AAC Rep1", ylabel="AAC Rep2")
    hexbin!(ax8, completePairs_qq.sorted_1, completePairs_qq.sorted_2, cellsize=cs, threshold=1, colormap=cl_map)
    ablines!(ax8, [0], [1])
    text!(ax8, posX, posY,text="rₛ=$(qqCorr_complete[1, :rₛ])\nr=$(qqCorr_complete[1, :r])", align=alignArg)

    ## Incomplete pairs
    incompletePairs_qq = filter([:viab_sd_1, :viab_sd_2] => (a, b) -> a < 20 && b < 20, paired_posterior)
    qqCorr_incomplete = correlationAnalysis(incompletePairs_qq[:,:sorted_1], incompletePairs_qq[:,:sorted_2])

    ax9 = Axis(f[3, 3], title="Posterior QQ\nIncomplete pairs(n=$(nrow(incompletePairs_median)))", xlabel="AAC Rep1", ylabel="AAC Rep2")
    hexbin!(ax9, incompletePairs_qq.sorted_1, incompletePairs_qq.sorted_2, cellsize=cs, threshold=1, colormap=cl_map)
    ablines!(ax9, [0], [1])
    text!(ax9, posX, posY,text="rₛ=$(qqCorr_incomplete[1, :rₛ])\nr=$(qqCorr_incomplete[1, :r])", align=alignArg)

    linkaxes!(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)
    return f
end

function plot_negativeAACposterior(negAAC_posterior::DataFrame, sample_expId::DataFrame)
    K = negAAC_posterior.exp_id |> unique |> length ## nb of experiments
    Kc = sum(sample_expId.group .== 1)
    Ki = sum(sample_expId.group .== 2)

    k = nrow(negAAC_posterior) ## nb of MCMC samples
    kc = sum(negAAC_posterior.viab_sd .≥ 20) ## nb of sample with -AAC but "complete curves"
    ki = sum(negAAC_posterior.viab_sd .< 20) ## nb of sample with -AAC and "incomplete curves"

    f = Figure(backgroundcolor="transparent", resolution = (500*3, 400*1));
    ax1 = Axis(f[1:2, 1], title="LDR vs. HDR - Samples with -AAC", xlabel="LDR", ylabel="HDR")
    ax2 = Axis(f[1:2, 2], title="AAC vs. SD - Samples with -AAC", xlabel="AAC", ylabel="Viab. SD (%)")
    ax3 = Axis(f[1, 3], title="Nb. of -AAC samples by exeriment - SD ≥ 20", ylabel="Nb. Exp.")
    ax4 = Axis(f[2, 3], title="Nb. of -AAC samples by exeriment - SD < 20", xlabel="Proportion MCMC samples", ylabel="Nb.Exp.")

    hexbin!(ax1, negAAC_posterior.LDR, negAAC_posterior.HDR, cellsize=0.5, threshold=1, colormap=["gray95", "gray35", "black"])
    ablines!(ax1, [0], [1])
    text!(ax1, 100, 80, text="Nb. Exp. = $K\nNb. MCMC samples = $k", align=[:center, :baseline])

    hexbin!(ax2, negAAC_posterior.aacRC, negAAC_posterior.viab_sd, cellsize=0.5, threshold=1, colormap=["gray95", "gray35", "black"])
    hlines!(ax2, [20])
    text!(ax2, -40, 30, text="\nNb. MCMC samples = $kc", align=[:center, :baseline])
    text!(ax2, -40, 10, text="\nNb. MCMC samples = $ki", align=[:center, :baseline])

    tmp = filter(:group => x -> x == 1, sample_expId)
    hist!(ax3, tmp.propIte, color=:black)
    text!(ax3, 0.2, 20, text="Nb. Exp. = $Kc", align=[:center, :baseline])
    tmp = filter(:group => x -> x == 2, sample_expId)
    hist!(ax4, tmp.propIte, bins=50, color=:black)
    text!(ax4, 0.2, 200, text="Nb. Exp. = $Ki", align=[:center, :baseline])

    hidexdecorations!(ax3, grid=false)
    linkxaxes!(ax3, ax4)

    return f
end