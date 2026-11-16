# Run one seasonal-forcing slice of the influenza parameter sweep.
# Called from Slurm; each array task writes into its own job-specific directory.
using CSV, DataFrames

coverage = parse(Int, get(ENV, "INFLUENZA_COVERAGE", "0"))
coverage in (30, 80) || error("INFLUENZA_COVERAGE must be 30 or 80")
coverage_campaign = coverage / 100

include(joinpath(@__DIR__, "..", "models", "influenza_model.jl"))

task_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
1 <= task_id <= length(forcing_vals) || error("Invalid Slurm array task ID: $(task_id)")
batch_id = get(ENV, "SLURM_ARRAY_JOB_ID", "")
occursin(r"^\d+$", batch_id) || error("Missing or invalid SLURM_ARRAY_JOB_ID")

forcing = forcing_vals[task_id]
forcing_tag = lpad(string(round(Int, 100 * forcing)), 3, '0')
part_dir = joinpath(output_dir, "forcing_parts", "coverage$(coverage)_batch$(batch_id)")
mkpath(part_dir)

manifest_file = joinpath(part_dir, "manifest_f$(forcing_tag).csv")
isfile(manifest_file) && error("Manifest exists: $(manifest_file)")

println("Coverage=$(coverage)%, forcing=$(forcing), task=$(task_id), batch=$(batch_id)")
df = run_influenza_sweep(;
    forcing_values=[forcing],
    output_directory=part_dir,
    output_suffix="f$(forcing_tag)"
)

expected_rows = length(beta_vals) * length(nu_vals) * length(first_age_range)
nrow(df) == expected_rows || error("Expected $(expected_rows) rows, got $(nrow(df))")

manifest = DataFrame(;
    batch_id=[batch_id],
    coverage=[coverage],
    force_frac=[forcing],
    n_beta=[length(beta_vals)],
    n_nu=[length(nu_vals)],
    n_ages=[length(first_age_range)],
    expected_rows=[expected_rows],
    actual_rows=[nrow(df)],
)
CSV.write(manifest_file, manifest)
println("Complete forcing slice: $(forcing); manifest: $(manifest_file)")
