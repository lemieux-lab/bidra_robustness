# BiDRA Robustness Demonstration
Code for the evaluation of bidra robustness and application to discrepancies analysis. We are using three large scale datasets (gCSI, CTRPv2 and GRAY) that are stored on lab servers.

## Do compounds characterization
To predict estimates with standard Marquardt-Levenberg: `julia curveFit.jl`. The inference is for the experiments of all three experiments. The results are stored on the server.

To infer posterior distribution with BiDRA: `./partionBiDRA.sh [nameOfDataset]`. Temporary diagnostics files are created for each batch. To merge them into one diagnotics file: `julia mergeDiagnostics.jl`. Posterior distributions, complete chains results and figures are stored on the server.


