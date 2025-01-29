# BiDRA Robustness Demonstration
Code for the evaluation of bidra robustness and application to discrepancies analysis.

## Datasets
Used datasets: Gray (public), gCSI (public), CTRPv2 (public) and IRIC (in-house).

The data from the three public datasets are access through PharmacoGX, using the `output_[nameOfDataset]_curves.ipynb`. Additional relevant files are also downloaded. An anonymized version of the IRIC dataset is made available.

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
```





To get correlation metrics (slope of linear fit, r2, Spearman and Pearson correlation coefficients) of viability responses (shared concentration): `julia viabilityCorrelation.jl [nameOfDataset]`. Correlation metrics are calculated for 3 within-dose replicates representation: all possible pairing, mean responses, and bootstraping of pairing. For the later, the mean, median and std of correlation metrics are returned. If there are no within-concentration replicates, only the first representation is used. Results are added to a CSV file.

To get correlation metrics (slope of linear fit, r2, Spearman and Pearson correlation coefficients) of metrics estimates (Marquardt-Levenberg): `julia MLcorrelation.jl [nameOfDataset]`. Results are added to a CSV file. 

To get correlation metrics of posterior distributions (BiDRA): `julia posteriorCorrelation.jl [nameOfDataset]`. Correlation metrics are calculated for 4 posterior representations: median, three-percentiles, random sampling and quantile-to-quantile. Results are added to various CSV files.