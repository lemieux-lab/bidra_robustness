# BiDRA Robustness Demonstration
Code for the evaluation of bidra robustness and application to discrepancies analysis.

## Datasets
Used datasets: Gray (public), gCSI (public), CTRPv2 (public) and IRIC (in-house).

The data from the three public datasets are access through PharmacoGX, using the `output_[nameOfDataset]_curves.ipynb`. Additional relevant files are also downloaded. An anonymized version of the IRIC dataset is made available. It is important to make sure that the first column of the `curve_info/` files is labelled `exp_id`.

As a first step, the data `CSV` of the public datasets are converted to `H5` files with `public_data/csvToH5.jl`. Experiments are also filtered for extreme values, as describe in the manuscript.  


## Robustness and Discrepancies analysis
Datasets used: Gray, gCSI and CTRPv2

### Compounds characterization
Posterior for each experiments are generated with `compound_characterization/bidra.jl`. The bash script `compound_characterization/partitionBiDRA.sh` allows to split a given dataset in N batch and run BiDRA simultaneously. Posterior are stored in the dataset `H5` files (`public_datasets/bidra/`). Batches' diagnostics files can me merged for a single dataset with `compound_characterization/mergeDiagnostics.jl`. The diagnostics and the batch timing statistics are stored in `_generated_data/`.

Once all three datasets have been imported and converted to `H5`, LM estimates for each experiments can calculated with `compound_characterization/curveFit.jl`. Results are saved in a single file in `public_datasets/all_julia_curveFit.csv`.

### Correlation analysis
1. Pairings of duplicated experiments are obtained with `correlation_metrics/generatesPairings.jl`. This code also produce a list of of pairings for more than two replicated experiments. Results are saved in `public_datasets/`.

2. Response consistency across biological replicates is assessed through calculation of correlation metrics for viability responses. Values are obtained with `correlation_metrics/viabilityCorrelation.jl` and saved in `_generated_data/viabCorrelations.csv`.

3. Efficiency metrics correlation is calculated with `correlation_metrics/posteriorCorrelation.jl` and `correlation_metrics/LMcorrelation.jl`. The results are respectively saved in `_generated_data/posteriorCorrelations.csv`, `_generated_data/medianCorrelations.csv` and `_generated_data/qqCorrelations.csv`, and in `_generated_data/mlCorrelations.csv`.

3. Correlations between randomly paired experiments are calculated with `correlation_metrics/runRandomPairings_BiDRA.jl` and `correlation_metrics/runRandomPairings_ML.jl`. The results are respectively saved in `_generated_data/bidraRandomCorrelation.csv` and in `_generated_data/mlRandomCorrelations.csv`.

### Figures
Figures are outputed in `_generated_figures` and illustrate results obtained from the compound characterization and the correlation analysis.

### Other analysis
Analysis results presented as supplementary materials can be generated with the scripts contained in `other_analysis/`. These analysis include the multi-replicates correlations (`multiRep_posteriorcorrelation.jl` and `multiRep_MLcorrelation.jl`), the across datasets correlations (`acrossDataset_singletonCorr.jl`), the AAC analysis (`aac_analysis.jl`) and the generation of the various plots (`generateFigure_supp.jl` and `aac_plot.jl`).

## SAR Analysis
The analysis of the IRIC dataset is stand-alone and is designed to be run within the `iric_dataset` directory. Posterior are generated with `bidra_inference.jl` and saved in `results`. The complete SAR analysis can be run with `sar_analysis.jl` and the resulting figures will be saved in `figures`.

## Arborescence

```
project
│   README.md
│   Manifest.toml    
│   Project.toml
|
└───_generated_data
|   └───tmp
|
└───_generated_figures
|   └───discrepancies_replicates
|   └───methods_comparison
|   └───models
|   └───robustness
|   └───viab_corr
|   └───supp_fig
|
└───public_datasets
|   |   csvToH5.jl
|   |   output_ctrpv2_curves.ipynb
|   |   output_gCSI_curves.ipynb
|   |   output_gray_curves.ipynb
|   |
|   └───bidra
|   └───cellAnnotations
|   └───curves_info
|   └───drugAnnotations
|
└───iric_dataset
|   |   IRIC_anonymized.csv
|   |   bidra_inference.jl
|   |   sar_analysis.jl
|   |   utils_iric.jl
|   |
|   └───results
|   └───data
```
