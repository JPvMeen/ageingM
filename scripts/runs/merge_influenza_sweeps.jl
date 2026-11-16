# Combine completed seasonal-forcing Slurm-array slices into two batch-specific CSVs.
# Usage: julia --project=. scripts/runs/merge_influenza_sweeps.jl 80 123456
using CSV, DataFrames

length(ARGS) == 2 || error("Usage: merge_influenza_sweeps.jl <coverage:30|80> <array_job_id>")
coverage = parse(Int, ARGS[1])
coverage in (30, 80) || error("Coverage must be 30 or 80")
batch_id = ARGS[2]
occursin(r"^\d+$", batch_id) || error("Array job ID must be numeric")

root = normpath(joinpath(@__DIR__, "..", ".."))
parts_dir = joinpath(root, "output", "forcing_parts", "coverage$(coverage)_batch$(batch_id)")
out_dir = joinpath(root, "output")
mkpath(out_dir)

forcing_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
forcing_tags = [lpad(string(round(Int, 100 * f)), 3, '0') for f in forcing_values]
prefix = "influenza_2025_v50_forcing_sweep_coverage$(coverage)"
key_cols = [:beta, :nu, :first_age, :force_frac]
expected_grids = nothing
validated_results = Dict{String,DataFrame}()
final_paths = Dict(outcome => joinpath(out_dir, "$(outcome)_$(prefix)_batch$(batch_id).csv") for outcome in ("mortality", "yll"))
for (outcome, path) in final_paths
    isfile(path) && error("Combined output already exists: $(path)")
end

for outcome in ("mortality", "yll")
    parts = DataFrame[]
    for (forcing, tag) in zip(forcing_values, forcing_tags)
        input_path = joinpath(parts_dir, "$(outcome)_$(prefix)_f$(tag).csv")
        manifest_path = joinpath(parts_dir, "manifest_f$(tag).csv")
        isfile(input_path) || error("Missing $(outcome) result: $(input_path)")
        isfile(manifest_path) || error("Missing completed-task manifest: $(manifest_path)")

        manifest = CSV.read(manifest_path, DataFrame; types=Dict(:batch_id=>String))
        nrow(manifest) == 1 || error("Invalid manifest: $(manifest_path)")
        meta = manifest[1, :]
        string(meta.batch_id) == batch_id || error("Batch ID mismatch in $(manifest_path)")
        meta.coverage == coverage || error("Coverage mismatch in $(manifest_path)")
        isapprox(meta.force_frac, forcing; atol=1e-12) || error("Forcing mismatch in $(manifest_path)")
        expected = meta.n_beta * meta.n_nu * meta.n_ages
        expected == meta.expected_rows == meta.actual_rows ||
            error("Incomplete task in $(manifest_path)")

        df = CSV.read(input_path, DataFrame)
        nrow(df) == expected || error("Unexpected row count in $(input_path)")
        all(isapprox.(df.force_frac, forcing; atol=1e-12)) ||
            error("Incorrect forcing values in $(input_path)")
        nrow(unique(df, key_cols)) == nrow(df) ||
            error("Duplicate parameter combinations in $(input_path)")

        # Each forcing amplitude must use the same full R0/nu/vaccination-age grid.
        grid = select(df, :beta, :nu, :first_age)
        sort!(grid, [:beta, :nu, :first_age])
        if expected_grids === nothing
            expected_grids = grid
        else
            isequal(grid, expected_grids) || error("Parameter grids differ in $(input_path)")
        end
        push!(parts, df)
    end

    combined = vcat(parts...; cols=:setequal)
    nrow(combined) == sum(nrow, parts) || error("Merge row count mismatch")
    nrow(unique(combined, key_cols)) == nrow(combined) ||
        error("Duplicate parameter combinations after merge")

    validated_results[outcome] = combined
end

# Both outcomes are fully validated before any final output is written.
for outcome in ("mortality", "yll")
    output_file = final_paths[outcome]
    combined = validated_results[outcome]
    CSV.write(output_file, combined)
    println("Saved $(output_file) ($(nrow(combined)) rows)")
end
