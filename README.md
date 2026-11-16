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

The code is organised into separate files for different simulation scenarios, allowing scenarios to be run and modified independently.

- `data/` — Input data, including the contact matrix.
- `scripts/` — Julia model and scenario files, as well as scripts for analysis and plotting.
- `output/` or `outputs/` — Generated simulation results, including mortality and YLL calculations.
- `figures/` — Generated plots and heatmaps.

### Scenarios

The simulations are organised into separate scenario files in `scripts/`, each focusing on different demographic assumptions or vaccination strategies.

**General scenarios**

- `ageing_general_2025.jl` — General SIRV model using mortality assumptions for 2025. Performs parameter sweeps over transmission rates, immunity waning rates, and vaccination start ages.
- `ageing_general_2050.jl` — Similar to the 2025 scenario, but using mortality assumptions for 2050 to investigate the effects of demographic changes on mortality and YLL.
- `ageing_general_CMfrequency.jl` — Scenario intended for investigating contact matrix frequency assumptions. Currently contains the same implementation as the 2050 scenario.

**Influenza scenarios**

- `ageing_influenza_coverage30.jl` — Influenza-like scenario with seasonal transmission, annual vaccination campaigns, and 30% vaccination coverage. Explores different transmission rates, immunity durations, and vaccination start ages.
- `ageing_influenza_coverage80.jl` — Same influenza-like scenario, but with 80% vaccination coverage, allowing comparison of the effects of higher vaccination coverage on mortality and YLL.

Each scenario can be run independently and produces scenario-specific output files.

## Running simulations

Simulations are run using Julia. Each scenario is defined in a separate file and can be executed independently.

From the repository root, run:

```bash
julia scripts/<scenario_file>.jl
```

Replace `<scenario_file>.jl` with the relevant scenario script.

The scenarios allow different parameter configurations and parameter sweeps to be run without modifying the other scenario files.

## Input and output

### Input

- Input data is the contact data in `data/contact_matrix_data_120.txt`. Notice that some other pathogen- and population-specific parameters are defined in the code (to be streamlined).
- Scenario-specific parameters are defined in separate Julia scenario files.

### Output

- Mortality and years-of-life-lost (YLL) results are saved as scenario-specific CSV files in the output directory.
- Various plots and time–age heatmaps are also produced and saved in the figures directory.
- Generated output files and figures are ignored by Git and can be regenerated by running the relevant scenarios.

## Dependencies

- DifferentialEquations, DiffEqCallbacks, LinearAlgebra, CSV, DataFrames, Plots, QuadGK, Distributions, NLsolve.

## Status

- This code is work in progress and comes with no warranty whatsoever.