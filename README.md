# Age-structured SIRV transmission model with population ageing

**Note:** This is work in progress and comes with no warranty!

This repository contains code for simulating an age-structured SIRV model with semi-realistic demographic projections for the Netherlands, Erlang-distributed immunity durations, age-targeted vaccination, seasonal forcing, and mortality/years-of-life-lost (YLL) calculations. The code is written in Julia and makes use of the *DifferentialEquations.jl* ecosystem. Related code under development based on the *MethodOfLines.jl* package is available on request. 

For enquiries, please contact:
- Jan Paul van Meenen: r.j.p.vanmeenen-2@umcutrecht.nl
- Michiel van Boven: r.m.vanboven-2@umcutrecht.nl
- Mirjam Kretzschmar: m.e.e.kretzschmar@umcutrecht.nl

## What does the program do?

- Sets up an age grid (0–120 years, default step 0.1y) and simulates dynamics over long horizons.
- Constructs a contact matrix using data in `data/contact_matrix_data_120.txt`. The contact matrix is scaled so that its largest eigenvalue equals 1. Hence R0 = beta/gamma. The code has been checked for threshold behaviour at R0 = 1 (using fixed demography, of course). The code also yields the correct final size when all contact matrix entries are identical.
- Implements Erlang-distributed immunity, both after vaccination and infection.
- Allows for season-to-season variation in waning rates by sampling from a (Gamma) distribution.
- Includes seasonal forcing of transmission, age-window vaccination campaigns, and age-dependent background mortality.
- Computes mortality and years of life lost from infections, using an infection fatality rate function and life expectancy integrals.
- Currently, it runs parameter sweeps over the transmission rate, the waning rates, and the vaccination start age. These are organised into separate scenario files. Notice that the code is not yet streamlined; care should be taken when adapting the code to your own needs.
- Saves mortality and years of life lost (YLL) and produces plots and time–age heatmaps. Notice that YLL is based on cohort life expectancy, NOT period life expectancy.

## Method

- The code discretises the age-structured partial differential equations (PDEs) using a first-order upwind scheme in age and standard ODE solvers in time.
- Vaccination is implemented as short campaigns at specific ages with a rate chosen to achieve a target coverage in a given age window (think: in three/six months, mimicking the duration of the campaign). Alternatives can easily be added.
- Mortality and YLL are obtained by combining age-specific infection deaths with conditional life expectancy.

## Directory structure

The repository separates shared model code, simulation runs, and figure scripts:

- `data/` — Input data, including the contact matrix.
- `scripts/models/` — Shared Julia equations and sweep logic: `ageing_general_model.jl` and `influenza_model.jl`.
- `scripts/runs/` — Short Julia scenario entry points for general ageing and influenza simulations.
- `scripts/figure/` — Figure-generation scripts named `Figure1.R` to `Figure7.R`, with superseded figure versions in `scripts/figure/archive/`.
- `output/` or `outputs/` — Generated simulation results, including mortality and YLL calculations.
- `figures/` — Generated plots and heatmaps.

### Scenarios

**General scenarios**

- `scripts/runs/ageing_general_2025.jl` — General SIRV model using mortality assumptions for 2025.
- `scripts/runs/ageing_general_2050.jl` — General SIRV model using mortality assumptions for 2050.
- `scripts/runs/ageing_general_CMfrequency.jl` — Currently identical in behaviour to the 2050 scenario.

These all use `scripts/models/ageing_general_model.jl`. The 2025 and 2050 run scripts provide distinct mortality offsets, R0 grids, and output names.

**Influenza scenarios**

- `scripts/runs/ageing_influenza_coverage30.jl` — Seasonal influenza parameter sweep with 30% vaccination coverage.
- `scripts/runs/ageing_influenza_coverage80.jl` — Seasonal influenza parameter sweep with 80% vaccination coverage.
- `scripts/runs/ageing_influenza_timeseries.jl` — Generates `time_series_<tag>_fva<age>_R0_<R0>_nu_<nu>.csv` and `time_age_<tag>_fva<age>_R0_<R0>_nu_<nu>.csv` using the same influenza model.

Both influenza coverage sweeps and the time-series run use `scripts/models/influenza_model.jl`.

## Running simulations

From the repository root, run the scenario of interest:

```bash
julia --project=. scripts/runs/ageing_general_2025.jl
julia --project=. scripts/runs/ageing_general_2050.jl
julia --project=. scripts/runs/ageing_influenza_coverage80.jl
julia --project=. scripts/runs/ageing_influenza_timeseries.jl
```

The `test_ageing_*.sh` Slurm scripts use the same `scripts/runs/` paths.

Figure generation is separate. Run R figure scripts from the repository root, so their `output/`, `figures/`, and `functions/` paths resolve correctly:

```bash
Rscript scripts/figure/Figure4.R
```

For the old general-scenario Julia diagnostic plots, an optional standalone `scripts/figure/ageing_general_diagnostics.jl` script is available. It is not executed by the general parameter sweeps.


### Parallel influenza forcing sweeps (Slurm)

For full 30% or 80% vaccination-coverage forcing sweeps, submit a Slurm array
from the repository root:

```bash
bash submit_influenza_array.sh 80
bash submit_influenza_array.sh 30
```

Each submission creates a six-task Slurm array (one task per seasonal forcing:
0, 0.05, 0.10, 0.15, 0.20, 0.30) with four Julia threads per task.
By default, at most three tasks are run concurrently. To lower concurrency:

```bash
INFLUENZA_MAX_PARALLEL=2 bash submit_influenza_array.sh 80
```

Each array job has its own batch ID, with intermediate CSVs and task manifests in
`output/forcing_parts/coverage80_batch<array_job_id>/` (or `coverage30`).
This prevents mixing different batches, and the merge script checks each
manifest, parameter grid, row count, forcing value, and duplicate keys.

A separate Slurm merge job is submitted with an `afterok` dependency:
**it only starts after every array task succeeds**. Its outputs are

```text
output/mortality_influenza_2025_v50_forcing_sweep_coverage80_batch<array_job_id>.csv
output/yll_influenza_2025_v50_forcing_sweep_coverage80_batch<array_job_id>.csv
```

The combined CSVs retain the existing model's column names and all six forcing
levels. The `batch<array_job_id>` suffix avoids replacing previous datasets.
Figure 5 and Figure 6 currently read older, fixed filenames; when using a new
combined dataset, update the R CSV path and filter to `force_frac == 0.05`
before computing relative mortality. Keep the no-vaccination baseline
(`first_age == 200`) within that forcing level.

The existing single-job scripts `test_ageing_influenza_coverage30.sh` and
`test_ageing_influenza_coverage80.sh` remain available as before. The Slurm
array and merge scripts do **not** generate figures or time-series CSVs.

### Parameterized influenza time-series simulations

Submit a single scenario with Slurm, specifying `R0`, vaccine-immunity
**waning rate** `nu` (per year), and the first vaccination age:

```bash
mkdir -p logs
sbatch test_ageing_influenza_timeseries.sh 1.5 0.5 91
```

Here `nu=0.5` means a mean vaccine-immunity duration of `1/nu = 2` years.
The optional fourth argument is an output tag (default `figure4_baseline`).
The default coverage is 80%, forcing is 0.05, and vaccination begins at year 200.

The example saves `output/time_series_figure4_baseline_fva91_R0_1p5_nu_0p5.csv`
and `output/time_age_figure4_baseline_fva91_R0_1p5_nu_0p5.csv`.
Numbers in file tags use `p` for decimal points. Each parameter combination
gets a distinct filename. Both Figure4.R and Figure4_diagnostic.R now include
the appropriate `nu` tag when reading CSV data.

The Julia script also accepts these parameters directly:

```bash
julia --project=. scripts/runs/ageing_influenza_timeseries.jl 1.5 0.5 91
```

## Input and output

### Input

- Contact data is in `data/contact_matrix_data_120.txt`.
- Scenario-specific parameters are in the run scripts, with shared equations and parameters in `scripts/models/`.

### Output

- Mortality, years-of-life-lost (YLL), and time-series results are saved as scenario-specific CSV files in `output/`.
- Figure-generation scripts in `scripts/figure/` write plots to `figures/`.
- Generated output files and figures are ignored by Git and can be regenerated independently.

## Dependencies

- DifferentialEquations, DiffEqCallbacks, LinearAlgebra, CSV, DataFrames, Plots, QuadGK, Distributions, NLsolve.

## Status

- This code is work in progress and comes with no warranty whatsoever.