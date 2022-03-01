include("utils.jl")

### Batch variables
dataset = ARGS[1]
Nbatch = parse(Int64, ARGS[2])

diagnostic_fn = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/"*datate*"_diagnostics.csv"
diagnostic_data = DataFrame(exp_id=[], batch=[], time=[], HDR=[], LDR=[], ic50=[], slope=[], σ=[])

for i in 1:Nbatch
    batch = i-1
    diagnosticTMP_fn = "/u/labellec/Desktop/bayesian_dose_response/bidra_robustness/_generated_data/TMP_diagnostics"*String(batch)*".csv"

    tmp = readCSV(diagnosticTMP_fn, false)
    append!(diagnostic_data)
end

CSV.write(diagnostic_fn, diagnostic_data, delim=",")