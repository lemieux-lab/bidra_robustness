using DataFrames, CSV, HDF5, Statistics, StatsBase

dts = DataFrame([
    (name="Gray", mle="gray", id=:exp_id),
    (name="gCSI", mle="gCSI", id=:expid), 
    (name="CTRPv2", mle="ctrpv2", id=:experimentIds)
])
dts.sym = Symbol.(dts.name)

function copy_and_sort(h, df, buf, m, i, j)
    copyto!(buf, h["/$(df[i,j])/chains"])
    sort!(buf, dims=1)
    m[i, :, :] .= Float32.(buf)
end

function do_all(mi, mj, metric_idx, spearman=false)
    res = zeros(eltype(mi), size(mi, 2))
    Threads.@threads for i ∈ axes(mi, 2)
        res[i] = r_swap_mc(mi[:, i, metric_idx], mj[:, i, metric_idx], 10_000, spearman)
    end
    return res
end

function prep_q_df(n)
    gen_qs(k) = round.(((1:k) .- 0.5) ./ k, digits=2)
    tmp_q = gen_qs(n)
    q_df = DataFrame(sym=Symbol.(tmp_q), q=tmp_q)
    return q_df
end

struct PairedMetrics{D <: AbstractDataFrame, M <: AbstractArray, B <: AbstractVector{Bool}}
    di::D
    dj::D
    mi::M
    mj::M
    chains_colName::Dict{Symbol, Int64}
    subs::Dict{Symbol, B}
end

function PairedMetrics(dt, randomize=false)
    h = h5open("public_datasets/bidra/$(dt.name)_complete.h5")
    
    df = identify_replicates(h, dt, randomize)
    
    is_complete = falses(size(df))
    for (i, row) in enumerate(eachrow(df))
        is_complete[i, :] .= (std(h["$(row.rep_1)/data"][:,2]) >= 20., std(h["$(row.rep_2)/data"][:,2]) >= 20.)
    end
    
    subs = Dict(
        :both_c => all(is_complete; dims=2) |> vec,
        :both_i => all(.~is_complete; dims=2) |> vec,
        :all => trues(nrow(df))
    )

    mle_data = CSV.read("public_datasets/all_julia_curveFit.csv", DataFrame; pool = true)

    nt(row) = (LDR=row.LDR, HDR=row.HDR, ic50=row.ic50, slope=row.slope, aac=row.aac, rmsd=row.rmsd, convergence=row.convergence)
    mle = Dict(Symbol(row.exp_id) => nt(row) for row ∈ eachrow(mle_data) if row.dataset == dt.mle)

    chains_colName = Dict(Symbol(s) => i for (i, s) ∈ enumerate(h["info/chains_colNames"][:]))

    di = DataFrame(mle[i] for i ∈ df[:,1])
    dj = DataFrame(mle[i] for i ∈ df[:,2])

    # To extract the chains...
    l_chain, l_param = size(h["/$(df[1,1])/chains"])

    mi = zeros(Float32, nrow(df), l_chain, l_param)
    mj = zeros(Float32, nrow(df), l_chain, l_param)
    buf = zeros(l_chain, l_param)

    for i ∈ 1:nrow(df)
        copy_and_sort(h, df, buf, mi, i, 1)
        copy_and_sort(h, df, buf, mj, i, 2)
    end

    return PairedMetrics(di, dj, mi, mj, chains_colName, subs)
end

function launch_analysis!(pm, dt, q_df, all_df, spearman=false)

    # Launch analysis

    metrics = [:LDR, :HDR, :ic50, :slope]
    sub_labels = [:all, :both_c, :both_i]

    res = Dict{Tuple{Symbol, Symbol}, Vector{Float32}}()
    res_lock = ReentrantLock()

    Threads.@threads for metric ∈ metrics
        for sub_label ∈ sub_labels
            tmp = do_all(pm.mi[pm.subs[sub_label],:,:], pm.mj[pm.subs[sub_label],:,:], pm.chains_colName[metric], spearman)
            lock(res_lock) do 
                res[(metric, sub_label)] = tmp
            end
        end
    end

    for metric ∈ metrics, sub_label ∈ sub_labels
        for row ∈ eachrow(q_df)
            row.sym == :med && continue
            push!(all_df, (
                dataset=Symbol(dt.name),
                method=:bidra_quantile,
                metric=metric,
                sub=sub_label,
                q=row.sym,
                r_swap=quantile(res[(metric, sub_label)], row.q)
            ))
        end

        push!(all_df, (
            dataset=Symbol(dt.name),
            method=:mle,
            metric=metric,
            sub=sub_label,
            q=Symbol("0.5"),
            r_swap=r_swap_mc(pm.di[pm.subs[sub_label], metric], pm.dj[pm.subs[sub_label], metric], 10_000, spearman)
        ))
    end
end

function prep_all_data(dts, n_quantile=15, randomize=false)
    q_df = prep_q_df(n_quantile)
    q_df.color = my_col.(q_df.q)

    dict_all_df = Dict{String, Any}()
    for cor_fn ∈ ["Pearson", "Spearman"]
        all_df = DataFrame()
        for dt ∈ eachrow(dts)
            pm = PairedMetrics(dt, randomize)
            launch_analysis!(pm, dt, q_df, all_df, cor_fn == "Spearman")
        end
        dict_all_df[cor_fn] = all_df
    end

    return dict_all_df, q_df
end