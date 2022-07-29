module Utils
export StrIndex, getindex, length
using HDF5, DataFrames

## Provide a 2-way indexing between string and int
struct StrIndex
    str2id::Dict{String, Int32}
    id2str::Vector{String}

    StrIndex(vs::Vector{String}) = new(Dict(vs[i] => i for i = 1:length(vs)), vs) # must be uniqued!
    StrIndex(ds::HDF5.Dataset) = StrIndex(ds[:])
end

## Indexing
Base.getindex(idx::StrIndex, s::String) = idx.str2id[s]
Base.getindex(idx::StrIndex, i::Integer) = idx.id2str[i]
Base.getindex(idx::StrIndex, v::AbstractVector{String}) = [idx[s] for s in v]
Base.getindex(idx::StrIndex, v::AbstractVector{<:Integer}) = [idx[i] for i in v]
Base.getindex(idx::StrIndex, df::AbstractDataFrame) = mapcols(col -> idx[col], df)
Base.length(idx::StrIndex) = length(idx.id2str)

## HDF5 IO
Base.setindex!(f::HDF5.File, s::StrIndex, k::String) = setindex!(f, s.id2str, k)

function Base.setindex!(f::HDF5.File, df::AbstractDataFrame, k::String)
    g = create_group(f, k)
    for (name, vec) in pairs(eachcol(df))
        g[String(name)] = vec
    end
end

function DataFrames.DataFrame(g::HDF5.Group)
    convert(p) = (p.first, p.second[:]) # To pull data from the HDF5 dataset
    return DataFrame(Dict(map(convert, pairs(g))))
end

end # module
