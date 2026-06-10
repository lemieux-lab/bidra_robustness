using DataFrames, HDF5
using CSV
using Statistics, StatsBase
using Distributions
using Gadfly, StatsPlots
using Cairo, Fontconfig
using CairoMakie

R = 10
data_prefix = "_generated_data_reviews";
synthetic_fn = data_prefix*"/synthetic_data_gaussian_noise_time_$(string(R))_repeats_oni.csv";

data = CSV.read(synthetic_fn, DataFrame);

# Remove the very first time measurements (repetition 1, sigma=1.0, subset_size=12 for both methods)
data = filter(row -> !(row.repetition == 1 && row.sigma == 1.0 && row.subset_size == 12), data)
println("Filtered out repetition 1 for σ=1.0, subset_size=12")

figure_prefix = "_generated_figures_reviews/synthetic_missing_doses/";
if !isdir(figure_prefix)
    mkpath(figure_prefix)
end

# Get unique sigma values
sigma_values = sort(unique(data.sigma))
println("Creating figure for sigma values: ", sigma_values)

# Create figure with 3 rows and 2 columns (one column per sigma)
fig = CairoMakie.Figure(backgroundcolor="transparent", size=(1200, 1000));

for (sigma_idx, sigma) in enumerate(sigma_values)
    # Filter data for this sigma
    sigma_data = filter(row -> row.sigma == sigma, data)
    bidra_sigma = filter(row -> row.method == "BiDRA", sigma_data)
    lsqfit_sigma = filter(row -> row.method == "LsqFit", sigma_data)
    
    # Get unique subset sizes and sort in descending order (12 to 2)
    subset_sizes = sort(unique(sigma_data.subset_size), rev=true)
    
    # Create categorical labels (strings) for x-axis
    subset_labels = string.(subset_sizes)
    
    # Column position for this sigma
    col = sigma_idx
    
    # Prepare data for boxplots and calculate medians
    bidra_medians = Float64[]
    lsqfit_medians = Float64[]
    ratio_medians = Float64[]
    
    # ===== Row 1: BiDRA Plot =====
    ax_bidra = CairoMakie.Axis(fig[1, col], 
        xlabel="", 
        ylabel="Time (seconds)", 
        title="σ=$(sigma) - BiDRA",
        xticks=(1:length(subset_sizes), subset_labels),
        xticklabelsvisible=false)
    
    for (idx, ss) in enumerate(subset_sizes)
        # Get all time values for this subset_size
        bidra_times = filter(row -> row.subset_size == ss && row.method == "BiDRA", sigma_data).time
        
        # Store medians for line plots
        push!(bidra_medians, median(bidra_times))
        
        # Create boxplot at categorical position
        CairoMakie.boxplot!(ax_bidra, fill(idx, length(bidra_times)), bidra_times, 
            width=0.5, color=(:blue, 0.5), whiskerwidth=0.5)
    end
    
    # Draw line connecting medians
    CairoMakie.lines!(ax_bidra, 1:length(subset_sizes), bidra_medians, 
        color=:blue, linewidth=2)
    
    # ===== Row 2: LsqFit Plot =====
    ax_lsqfit = CairoMakie.Axis(fig[2, col], 
        xlabel="", 
        ylabel="Time (seconds)", 
        title="σ=$(sigma) - LsqFit",
        xticks=(1:length(subset_sizes), subset_labels),
        xticklabelsvisible=false)
    
    for (idx, ss) in enumerate(subset_sizes)
        # Get all time values for this subset_size
        lsqfit_times = filter(row -> row.subset_size == ss && row.method == "LsqFit", sigma_data).time
        
        # Store medians for line plots
        push!(lsqfit_medians, median(lsqfit_times))
        
        # Create boxplot at categorical position
        CairoMakie.boxplot!(ax_lsqfit, fill(idx, length(lsqfit_times)), lsqfit_times, 
            width=0.5, color=(:red, 0.5), whiskerwidth=0.5)
    end
    
    # Draw line connecting medians
    CairoMakie.lines!(ax_lsqfit, 1:length(subset_sizes), lsqfit_medians, 
        color=:red, linewidth=2)
    
    # ===== Row 3: Ratio Plot (with log scale) =====
    ax_ratio = CairoMakie.Axis(fig[3, col], 
        xlabel="Subset Size (number of doses)", 
        ylabel="BiDRA/LsqFit", 
        title="σ=$(sigma) - Time Ratio",
        xticks=(1:length(subset_sizes), subset_labels),
        yscale=log10)
    
    for (idx, ss) in enumerate(subset_sizes)
        # Get BiDRA and LsqFit data for this subset, sorted by repetition
        bidra_subset = filter(row -> row.subset_size == ss && row.method == "BiDRA", sigma_data)
        lsqfit_subset = filter(row -> row.subset_size == ss && row.method == "LsqFit", sigma_data)
        
        # Sort by repetition to ensure proper matching
        sort!(bidra_subset, :repetition)
        sort!(lsqfit_subset, :repetition)
        
        # Calculate ratio for each repetition
        ratios = bidra_subset.time ./ lsqfit_subset.time
        
        # Store medians for line plots
        push!(ratio_medians, median(ratios))
        
        # Create boxplot at categorical position
        CairoMakie.boxplot!(ax_ratio, fill(idx, length(ratios)), ratios, 
            width=0.5, color=(:purple, 0.5), whiskerwidth=0.5)
    end
    
    # Draw line connecting medians
    CairoMakie.lines!(ax_ratio, 1:length(subset_sizes), ratio_medians, 
        color=:purple, linewidth=2)
    
    # Add horizontal line at 1.0 (equal times)
    CairoMakie.hlines!(ax_ratio, [1.0], color=:black, linewidth=1, linestyle=:dash)
    
    # Link x-axes of all three plots in this column
    CairoMakie.linkxaxes!(ax_bidra, ax_lsqfit, ax_ratio)
end

# Save the figure
output_fn = figure_prefix*"timing_comparison_$(string(R))_repeats_oni.pdf"
CairoMakie.save(output_fn, fig)

println("Figure saved to: ", output_fn)

# Print summary statistics for each sigma
println("\n=== Summary Statistics ===")
for sigma in sigma_values
    sigma_data = filter(row -> row.sigma == sigma, data)
    bidra_sigma = filter(row -> row.method == "BiDRA", sigma_data)
    lsqfit_sigma = filter(row -> row.method == "LsqFit", sigma_data)
    
    println("\nσ = $(sigma):")
    println("  BiDRA  - Mean: $(round(mean(bidra_sigma.time), digits=3))s, Median: $(round(median(bidra_sigma.time), digits=3))s")
    println("  LsqFit - Mean: $(round(mean(lsqfit_sigma.time), digits=6))s, Median: $(round(median(lsqfit_sigma.time), digits=6))s")
    println("  Speedup (LsqFit vs BiDRA): $(round(mean(bidra_sigma.time) / mean(lsqfit_sigma.time), digits=1))x")
end
