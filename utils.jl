using DataFrames
using CSV

function readCSV(fn) 
    csv_file = CSV.File("$fn")
    return DataFrame(csv_file); 
end

function replaceChar(df, col)
    df[!, col] .= replace.(df[:, col], " " => "_")
    df[!, col] .= replace.(df[:, col], ":" => "_")
    df[!, col] .= replace.(df[:, col], "/" => "_")
    return df
end

