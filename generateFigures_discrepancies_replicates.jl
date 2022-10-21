push!(LOAD_PATH, "Utils/")
using DataFrames
using Distributions, Statistics, StatsBase
using Gadfly, StatsPlots
using Cairo, Fontconfig
using Utils

### Global Var
figure_prefix = "_generated_figures/discrepancies_replicates/"

function correlationBarPlot(df::DataFrame, t::String, coef::Symbol)
    p=Gadfly.plot(df, x=:description, y=coef, color=:dataset, xgroup=:param, ygroup=:dataset,
                  Geom.subplot_grid( Geom.bar(position=:dodge)),
                  Guide.title(t))
    display(p)
end

function addMethods(df::DataFrame, method::String)
    df[!,:method] = repeat([method], nrow(df))
    return df
end

### Import correlation
ml_correlation = readCSV("_generated_data/mlCorrelations.csv", true)
median_correlation = readCSV("_generated_data/medianCorrelations.csv", true)
ci_correlation = readCSV("_generated_data/ciCorrelations.csv", true,)
posterior_correlation = readCSV("_generated_data/posteriorCorrelations.csv", true)
qq_correlation = readCSV("_generated_data/qqCorrelations.csv", true, )

### Plot correlations of individual methods
Gadfly.set_default_plot_size(5inch, 4inch)
correlationBarPlot(ml_correlation, "Marqaurdt-Levenberg Correlations", :rₛ)
correlationBarPlot(median_correlation, "BiDRA's median Correlations", :rₛ)
correlationBarPlot(qq_correlation, "BiDRA's QQ Correlations", :rₛ)

### Combine all correlations methods results and plot
ml_correlation = addMethods(ml_correlation, "ML")
median_correlation = addMethods(median_correlation, "Median")
qq_correlation = addMethods(qq_correlation, "QQ")
all_correlation = vcat(ml_correlation, median_correlation, qq_correlation)

### Plot all pairing results
Gadfly.set_default_plot_size(6inch, 4inch)
all_correlation_pairs = filter([:description, :param] => (x, y) -> x ∈ ["all pairs", "converged pairs"] && y != "aac", all_correlation)
p = Gadfly.plot(all_correlation_pairs, x=:method, y=:rₛ, color=:description, xgroup=:param, ygroup=:dataset, 
            Geom.subplot_grid(Geom.bar(position=:dodge)),
            Guide.title("All Pairs correlations"))
draw(PDF("$figure_prefix"*"all_dt_correlations_spearman.pdf", 6inch, 4inch), p)

p = Gadfly.plot(all_correlation_pairs, x=:method, y=:r, color=:description, xgroup=:param, ygroup=:dataset, 
            Geom.subplot_grid(Geom.bar(position=:dodge)),
            Guide.title("All Pairs correlations"))
draw(PDF("$figure_prefix"*"all_dt_correlations_pearson.pdf", 6inch, 4inch), p)

### Plot correlation by "completeness"
corr_subset = filter(:description => x -> x ∉ ["mixte pairs", "converged pairs"], all_correlation)
dt = "gCSI"
corr_subset_dt = filter(:dataset => x -> x == dt, corr_subset)

Gadfly.set_default_plot_size(10inch, 8inch)
p = Gadfly.plot(corr_subset, x=:description, y=:rₛ, color=:method, ygroup=:param, xgroup=:dataset, Geom.subplot_grid(Geom.bar(position=:dodge)))
draw(PDF("$figure_prefix"*"all_dt_correlations_grouped_spearman.pdf", 10inch, 8inch), p)

Gadfly.set_default_plot_size(10inch, 8inch)
p = Gadfly.plot(corr_subset, x=:description, y=:r, color=:method, ygroup=:param, xgroup=:dataset, Geom.subplot_grid(Geom.bar(position=:dodge)))
draw(PDF("$figure_prefix"*"all_dt_correlations_grouped_pearson.pdf", 10inch, 8inch), p)

p = Gadfly.plot(corr_subset, x=:r, y=:rₛ, ygroup=:dataset, xgroup=:param, color=:description, shape=:method,
                    Geom.subplot_grid(Geom.point, Geom.abline()))
draw(PDF("$figure_prefix"*"all_dt_correlations_spearman_pearson.pdf", 10inch, 8inch), p)


### Correlations of random pairings
bidra_randomRep = readCSV("_generated_data/bidraRandomCorrelation.csv", true)
bidra_randomRep_qq = filter(:method => x -> x == "qq", bidra_randomRep)
bidra_randomRep_qq = filter(:param => x -> x != "aac", bidra_randomRep_qq)

ml_randomRep = readCSV("_generated_data/mlRandomCorrelations.csv", true)
ml_randomRep[!, :description] = ml_randomRep[:, :method]
ml_randomRep[!, :method] = repeat(["ML"], nrow(ml_randomRep))

randomRep_compare = vcat(bidra_randomRep_qq, ml_randomRep)

Gadfly.set_default_plot_size(10inch, 6inch)
p = Gadfly.plot(randomRep_compare, x=:description, y=:rₛ , color=:dataset, xgroup=:param, ygroup=:method, 
            Geom.subplot_grid(Geom.boxplot()),
            Guide.title("Random Pairings All dataset"))
#display(p)
draw(PDF("$figure_prefix"*"random_pairings_spearman.pdf", 10inch, 6inch), p)

p = Gadfly.plot(randomRep_compare, x=:description, y=:r , color=:dataset, xgroup=:param, ygroup=:method, 
            Geom.subplot_grid(Geom.boxplot()),
            Guide.title("Random Pairings All dataset"))
#display(p)
draw(PDF("$figure_prefix"*"random_pairings_pearson.pdf", 10inch, 6inch), p)

### Plot pairing correlation by param x dataset x method x description
dataset = ["gray", "gCSI", "ctrpv2"]
bidra_params = [:LDR, :HDR, :ic50, :slope, :aac]
metrics_bounds = Dict(:HDR=>[-50,150], :LDR=>[70,150], :ic50=>[-10,10], :slope=>[0,10], :aac=>[0,100])
scaleColor = Scale.lab_gradient("gray95","black")

for i in 1:length(dataset)
    dt = dataset[i]
    println("### $dt")

    println("1. Get data and metrics")
    print("--> List of expId ")
    @time expId_list = getExpId_h5(dt);

    print("--> StrIndex ")
    @time si = Utils.StrIndex(expId_list);

    println("--> Raw data ")
    @time data_df = getRawData_h5(dt, false, si);

    println("--> Calculate Viab. SD")
    sd_df = combine(groupby(data_df, :exp_id), :Viability => std => :std_viability)
    expId_complete = sd_df[sd_df.std_viability .>= 20, :exp_id]
    expId_incomplete = sd_df[sd_df.std_viability .< 20, :exp_id]

    println("2. Get pairings")
    println("---> all")
    pairings_df = getPairings_h5(dt, si)
    println("---> complete a.k.a both exp. with SD ≥ 20")
    pairingComplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_complete && y ∈ expId_complete, pairings_df)
    println("---> incomplete a.k.a both exp with SD < 20")
    pairingIncomplete_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete && y ∈ expId_incomplete, pairings_df)
    println("---> mixte a.k.a at least one exp. has SD < 20")
    pairingMixte_df = filter([:rep_1, :rep_2] => (x, y) -> x ∈ expId_incomplete || y ∈ expId_incomplete, pairings_df)

    println("---> converge a.k.a. both exp. converge for LM")
    mlPaired_df = getMLestimates(dt, si, pairings_df)
    pairingConverge_df = filter([:convergence_rep1, :convergence_rep2] => (a, b) -> a + b == 2, mlPaired_df)

    println("3. Pairings descriptors")
    pairs_desc = [pairings_df, pairingComplete_df, pairingIncomplete_df, pairingMixte_df, pairingConverge_df]
    pairs_names = ["all pairs", "complete pairs", "incomplete pairs", "mixte pairs", "converged pairs"]
    N = [nrow(df) for df in pairs_desc]
    
    println("4. Get posterior pairings")
    all_type_rep = DataFrame()
    for j in 1:length(N)
        println("---> $j/",length(N))
        tmp1 = getPosterior_h5(dt, false, si, Array(pairs_desc[j].rep_1))
        rename!(tmp1, map(x -> "$x"*"_rep1", names(tmp1)))
        tmp2 = getPosterior_h5(dt, false, si, Array(pairs_desc[j].rep_2))
        rename!(tmp2, map(x -> "$x"*"_rep2", names(tmp2)))
        tmp = hcat(tmp1, tmp2)

        tmp[!, :description] = repeat([pairs_names[j]], nrow(tmp))

        all_type_rep = vcat(all_type_rep, tmp)
    end

    println("5. Plot posterior pairs by metrics")
    Gadfly.set_default_plot_size(10inch, 4inch)
    for pr in bidra_params
        println(pr)
        pr_rep1 = Symbol(string(pr)*"_rep1")
        pr_rep2 = Symbol(string(pr)*"_rep2")

        prQQ = combine(groupby(all_type_rep, [:exp_id_rep1, :description]), pr_rep1 => sort => :sorted_1, pr_rep2 => sort => :sorted_2)
        prQQ_zoom = filter([:sorted_1, :sorted_2] => (a, b) -> metrics_bounds[pr][1] ≤ a ≤ metrics_bounds[pr][2] && metrics_bounds[pr][1] ≤ b ≤ metrics_bounds[pr][2], prQQ)

        println("---> After zooming")
        pairs_after_zoom = combine(groupby(prQQ_zoom, :description), :exp_id_rep1 => length∘unique => :n)

        p = Gadfly.plot(prQQ_zoom, x=:sorted_1, y=:sorted_2,
                        xgroup=:description, 
                        Geom.subplot_grid(Geom.hexbin()),
                        Scale.color_continuous(colormap=scaleColor, minvalue=1),
                        Guide.title(string("$dt $pr pairing posterior (N =", N, ") (n=,", Array(pairs_after_zoom.n))));
        #display(p)
        draw(PDF("$figure_prefix"*"$dt"*"_"*string(pr)*"_posterior.pdf", 10inch, 4inch), p)
    end

    println("6. Get LM pairings")
    all_type_rep_LM = DataFrame()
    for j in 1:length(N)
        println("---> $j/",length(N))
        tmp = getMLestimates([dt], pairs_desc[j][:, [:rep_1, :rep_2]])
        tmp[!, :description] = repeat([pairs_names[j]], nrow(tmp))
        all_type_rep_LM = vcat(all_type_rep_LM, tmp)
    end

    println("7. Plot LM pairs by metrics")
    Gadfly.set_default_plot_size(10inch, 4inch)
    for pr in bidra_params
        println(pr)
        pr_rep1 = Symbol(string(pr)*"_rep1")
        pr_rep2 = Symbol(string(pr)*"_rep2")

        prLM_zoom = filter([pr_rep1, pr_rep2] => (a, b) -> metrics_bounds[pr][1] ≤ a ≤ metrics_bounds[pr][2] && metrics_bounds[pr][1] ≤ b ≤ metrics_bounds[pr][2], all_type_rep_LM)

        println("---> After zooming")
        pairs_after_zoom = combine(groupby(prLM_zoom, :description), :rep_1 => length∘unique => :n)

        p = Gadfly.plot(prLM_zoom, x=pr_rep1, y=pr_rep2,
                        xgroup=:description, 
                        Geom.subplot_grid(Geom.hexbin()),
                        Scale.color_continuous(colormap=scaleColor, minvalue=1),
                        Guide.title(string("$dt $pr pairing posterior (N =", N, ") (n=,", Array(pairs_after_zoom.n))));
        #display(p)
        draw(PDF("$figure_prefix"*"$dt"*"_"*string(pr)*"_LM.pdf", 10inch, 4inch), p)
    end
end