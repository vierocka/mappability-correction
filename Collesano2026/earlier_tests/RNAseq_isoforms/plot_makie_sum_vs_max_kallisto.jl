using CSV, DataFrames, CairoMakie

datasets = Dict(
    "HCC1937" => "HCC1937/kalisto",
    "JY"      => "JY",
    "HCT116"  => "HCT116"
)

sumfile = "peptide_9mers_ExprSum.tsv"
maxfile = "peptides_9mer_max.tsv"

function read_merge(folder)
    sumdf = CSV.read(joinpath(folder, sumfile), DataFrame; delim='\t')
    maxdf = CSV.read(joinpath(folder, maxfile), DataFrame; delim='\t')

    # col1 = kmer, col3 = expr
    sumdf = select(sumdf, 1 => :kmer, 3 => :sum_expr)
    maxdf = select(maxdf, 1 => :kmer, 3 => :max_expr)

    df = innerjoin(sumdf, maxdf, on=:kmer)

    x = log10.(Float64.(df.sum_expr) .+ 1.0)
    y = log10.(Float64.(df.max_expr) .+ 1.0)

    return x, y
end

mkpath("plots")

for (name, folder) in datasets
    x, y = read_merge(folder)

    fig = Figure(size = (650, 650))
    ax = Axis(fig[1,1], title=name,
              xlabel="log10(Sum est_counts + 1)",
              ylabel="log10(Max est_counts + 1)")

    # For large N, use density via hexbin-style: Makie has `hexbin` in MakieExtensions,
    # but a simple scatter may be heavy. We'll downsample if huge.
    n = length(x)
    if n > 500_000
        idx = rand(1:n, 500_000)
        scatter!(ax, x[idx], y[idx], markersize=1)
    else
        scatter!(ax, x, y, markersize=1)
    end

    lo = max(minimum(x), minimum(y))
    hi = min(maximum(x), maximum(y))
    lines!(ax, [lo, hi], [lo, hi], linestyle=:dash)

    save("plots/$(name)_sum_vs_max_makie.png", fig)
end
