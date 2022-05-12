# BiDRA Robustness Demonstration
Code for the evaluation of bidra robustness and application to discrepancies analysis. We are using three large scale datasets (gCSI, CTRPv2 and GRAY) that are stored on lab servers.

## Do compounds characterization
To predict estimates with standard Marquardt-Levenberg: `julia curveFit.jl`. The inference is for the experiments of all three experiments. The results are stored on the server.

To infer posterior distribution with BiDRA: `./partionBiDRA.sh [nameOfDataset]`. Temporary diagnostics files are created for each batch. To merge them into one diagnotics file: `julia mergeDiagnostics.jl`. Posterior distributions, complete chains results and figures are stored on the server.


## Do correlation analysis
Pairings of duplicated experiments are listed in `[nameOfDataset]_rep2_pairing.csv` files. To get those lists: `julia generatesPairings.jl [nameOfDataset]`.

To get correlation metrics (slope of linear fit, r2, Spearman and Pearson correlation coefficients) of paramters estimates (Marquardt-Levenberg): `julia MLcorrelation.jl [nameOfDataset]`. Results are added to a CSV file. 

To get correlation metrics of posterior distributions (BiDRA): `julia posteriorCorrelation.jl [nameOfDataset]`. Correlation metrics are calculated for 4 posterior representations: median, three-percentiles, random sampling and quantile-to-quantile. Results are added to various CSV files.