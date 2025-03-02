include("utils.jl")

### Batch variables
dataset = ARGS[1]
Nbatch = parse(Int64, ARGS[2])

diagnostic_fn = "_generated_data/"*dataset*"_diagnostics.csv"
diagnostic_data = DataFrame(exp_id=[], batch=[], time=[], HDR=[], LDR=[], ic50=[], slope=[], σ=[])

for i in 1:Nbatch
    batch = i-1
    diagnosticTMP_fn = "_generated_data/tmp/TMP_diagnostics"*string(batch)*".csv"
    
    tmp = readCSV(diagnosticTMP_fn, false)
    rename!(tmp, [:exp_id, :batch, :time, :HDR, :LDR, :ic50, :slope, :σ])
    append!(diagnostic_data, tmp)
end

CSV.write(diagnostic_fn, diagnostic_data, delim=",")