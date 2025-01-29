# BiDRA Robustness Demonstration
Code for the evaluation of bidra robustness and application to discrepancies analysis.

## Access data
Used datasets: Gray (public), gCSI (public), CTRPv2 (public) and IRIC (in-house).

The data from the three public datasets are access through PharmacoGX, using the `output_[nameOfDataset]_curves.ipynb`. Additional relevant files are also downloaded. An anonymized version of the IRIC dataset is made available.

## Do compounds characterization
To predict estimates with standard Marquardt-Levenberg: `julia curveFit.jl`. The inference is for the experiments of all three experiments. The results are stored on the server.

To infer posterior distribution with BiDRA: `./partionBiDRA.sh [nameOfDataset]`. Temporary diagnostics files are created for each batch. To merge them into one diagnotics file: `julia mergeDiagnostics.jl`. Posterior distributions, complete chains results and figures are stored on the server.


## Do correlation analysis
Pairings of duplicated experiments are listed in `[nameOfDataset]_rep2_pairing.csv` files. To get those lists: `julia generatesPairings.jl [nameOfDataset]`.

To get correlation metrics (slope of linear fit, r2, Spearman and Pearson correlation coefficients) of viability responses (shared concentration): `julia viabilityCorrelation.jl [nameOfDataset]`. Correlation metrics are calculated for 3 within-dose replicates representation: all possible pairing, mean responses, and bootstraping of pairing. For the later, the mean, median and std of correlation metrics are returned. If there are no within-concentration replicates, only the first representation is used. Results are added to a CSV file.

To get correlation metrics (slope of linear fit, r2, Spearman and Pearson correlation coefficients) of metrics estimates (Marquardt-Levenberg): `julia MLcorrelation.jl [nameOfDataset]`. Results are added to a CSV file. 

To get correlation metrics of posterior distributions (BiDRA): `julia posteriorCorrelation.jl [nameOfDataset]`. Correlation metrics are calculated for 4 posterior representations: median, three-percentiles, random sampling and quantile-to-quantile. Results are added to various CSV files.