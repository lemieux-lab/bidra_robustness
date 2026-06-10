using DataFrames, HDF5
using CSV
using Statistics, StatsBase
using Distributions
using Gadfly, StatsPlots
using Cairo, Fontconfig
using CairoMakie

include("utils.jl")

data_prefix = "_generated_data_reviews/correlation_reorder_pairings"
bidra_params = ["LDR", "HDR", "ic50", "slope"];
method_color = Dict("ML" => :black, "QQ" => :blue);

figure_prefix = "_generated_figures_reviews/pairs_ordering_correlation/";
if !isdir(figure_prefix)
    mkpath(figure_prefix)
end

### New random pairings correlation results
## Import correlation
ml_correlation = readCSV("$(data_prefix)/mlCorrelations.csv", true);
qq_correlation = readCSV("$(data_prefix)/qqCorrelations.csv", true);

## Add methods column
ml_correlation[:, :method] = repeat(["ML"], nrow(ml_correlation));
qq_correlation[:, :method] = repeat(["QQ"], nrow(qq_correlation));

## Combine all correlations methods results and plot
all_correlation = vcat(ml_correlation, qq_correlation);
unique_description = unique(all_correlation.description)
unique_correlation = unique(all_correlation.method)

all_correlation[:, :description_int] = [findfirst(unique_description .== dt) for dt in all_correlation.description];
all_correlation[:, :method_int] = [findfirst(unique_correlation .== m) for m in all_correlation.method];
all_correlation[:, :method_color] = [method_color[m] for m in all_correlation.method];

### Original correlations
original_ml_correlation = readCSV("_generated_data/mlCorrelations.csv", true)
original_qq_correlation = readCSV("_generated_data/qqCorrelations.csv", true)
original_ml_correlation[:, :method] = repeat(["ML"], nrow(original_ml_correlation));
original_qq_correlation[:, :method] = repeat(["QQ"], nrow(original_qq_correlation));

original_all_correlation = vcat(original_ml_correlation, original_qq_correlation);
original_all_correlation[:, :description_int] = [findfirst(unique_description .== dt) for dt in original_all_correlation.description];
original_all_correlation[:, :method_int] = [findfirst(unique_correlation .== m) for m in original_all_correlation.method];
original_all_correlation[:, :method_color] = [method_color[m] for m in original_all_correlation.method];

datasets = unique(all_correlation.dataset);
param_names = unique(all_correlation.param);

corr_metric = "rₛ";
for corr_metric in ["r", "rₛ"]
    ## Plot all pairing results
    fig = CairoMakie.Figure(backgroundcolor="white", size=(200*3, 200*4));

    for dt in datasets
        col = findfirst(datasets .== dt)
        println("Processing dataset: $dt")

        for pr in param_names
            row = findfirst(param_names .== pr)
            sub_df = filter(row -> row.dataset == dt && row.param == pr, all_correlation)
            sub_df_original = filter(row -> row.dataset == dt && row.param == pr && row.description_int !== nothing, original_all_correlation)

            x_val = sub_df.description_int
            x_val_original = Int.(sub_df_original.description_int)
            y_val = sub_df[:, Symbol(corr_metric)]
            y_val_original = sub_df_original[:, Symbol(corr_metric)]

            println("  Processing parameter: $pr")
            summary_stats = combine(
                groupby(sub_df, [:dataset, :description, :method]),
                :rₛ => mean => :rₛ_mean,
                :rₛ => std => :rₛ_std
            )
            println("    Summary statistics:")
            for row in eachrow(summary_stats)
                println("      $(row.method) ($(row.description)): $(row.rₛ_mean) ± $(row.rₛ_std)")
            end

            ax = CairoMakie.Axis(fig[row, col], title="$pr ($dt)", xlabel="Pairings", ylabel="$corr_metric",)
            CairoMakie.barplot!(ax, x_val_original, y_val_original, dodge=sub_df_original.method_int, color=sub_df_original.method_color)
            CairoMakie.boxplot!(ax, x_val, y_val, dodge=sub_df.method_int, color=sub_df.method_color)
            CairoMakie.ylims!(ax, 0.0, 1.0)
            CairoMakie.hlines!(ax, [0.5, 0.75, 1.0], color=[:red, :yellow, :green], linestyle=:dash)
        end
    end

    display(fig)

    fn = "$(figure_prefix)/$(corr_metric)_correlations_with_random_pair_orders.pdf"
    CairoMakie.save(fn, fig)
end