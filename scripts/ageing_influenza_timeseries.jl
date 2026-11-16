# One influenza simulation for Figure 4 (independent of the parameter sweeps).
# Change parameters and output_tag to generate a different scenario.
coverage_campaign = 0.80
include(joinpath(@__DIR__, "influenza_model.jl"))

R0_single       = 3.0
forcing_single  = 0.05
nu_duration     = 2.0       # mean vaccine immunity duration in years
first_age_s     = 60.0
t_vac_start     = 200.0
output_tag      = "figure4_baseline"

β_single = gamma * R0_single
ν_single = 1.0 / nu_duration

mkpath(output_dir)
mkpath(figures_dir)

mask_s = age_mask_from_first(age_points, first_age_s, end_age)
vac_s  = VacSched(mask_s, rate_in_window, camp_width, camp_phase, t_vac_start)

# solve single scenario with finer output for export: 20 time points per year
cache_s = ODECache(n_age, age_points, da, N0, gamma, forcing_single)
dt_out  = 1.0 / 20.0             # 20 points / year
prob_s  = ODEProblem(f_pde!, u0, (tmin, tmax),
                     (β_single, ν_single, contact_matrix, cache_s, vac_s, ALPHA_SCHED))
sol_s   = solve(prob_s, Tsit5(); reltol=1e-4, abstol=1e-6,
                callback=PositiveDomain(), saveat=dt_out, dense=false)

# unpack solution matrices
t     = sol_s.t
all_u = hcat(sol_s.u...)
S_mat  = all_u[1:n_age,                 :]'
I_mat  = all_u[n_age+1:2*n_age,         :]'
R1_mat = all_u[2*n_age+1:3*n_age,       :]'
V1_mat = all_u[3*n_age+1:4*n_age,       :]'
R2_mat = all_u[4*n_age+1:5*n_age,       :]'
V2_mat = all_u[5*n_age+1:6*n_age,       :]'
M_mat  = gamma .* (I_mat .* (ifr.(age_points))')  # instantaneous death *rate* per age bin

# totals over time (aggregated over age, using da as age-bin width)
tot_S = vec(sum(S_mat,  dims=2)) .* da
tot_I = vec(sum(I_mat,  dims=2)) .* da
tot_R = vec(sum(R1_mat .+ R2_mat, dims=2)) .* da
tot_V = vec(sum(V1_mat .+ V2_mat, dims=2)) .* da
tot_N = tot_S .+ tot_I .+ tot_R .+ tot_V

# total mortality incidence over time (aggregated over age)
mort_inc_total = vec(sum(M_mat, dims=2)) .* da   # integral over age

# export last 100 years
t_start_export = tmax - 100.0           # 100 years before tmax
idx_sel = findall(t .>= t_start_export)

# with dt_out=1/20, this gives 2001 points
if length(idx_sel) > 100 * 20
    idx_sel = idx_sel[2:end]
end

t_sel    = t[idx_sel]
S_sel    = tot_S[idx_sel]
I_sel    = tot_I[idx_sel]
R_sel    = tot_R[idx_sel]
V_sel    = tot_V[idx_sel]
N_sel    = tot_N[idx_sel]
M_sel    = mort_inc_total[idx_sel]

df_time = DataFrame(
    t     = t_sel,
    N     = N_sel,
    S     = S_sel,
    I     = I_sel,
    R     = R_sel,
    V     = V_sel,
    M_inc = M_sel,
)

CSV.write(joinpath(output_dir,
          "time_series_$(output_tag).csv"),
          df_time)

# take same time points (last 100 years, 20/year)
I_sel_full = I_mat[idx_sel, :]

# one age value per year: 120 ages (0..119 years)
bins_per_year = Int(round(1 / da))   # e.g. 10 if da = 0.1
age_inds_year = 1:bins_per_year:n_age
ages_year     = age_points[age_inds_year]   # ~ 0,1,2,...,119

I_ta = I_sel_full[:, age_inds_year]  # (time × 120)

df_I = DataFrame()
df_I.t = t_sel
for (j, age) in enumerate(ages_year)
    colname = Symbol("age_$(Int(round(age)))")
    df_I[!, colname] = I_ta[:, j]
end

CSV.write(joinpath(output_dir,
          "time_age_$(output_tag).csv"),
          df_I)

plt1 = plot(t, tot_S, label="S")
plot!(plt1, t, tot_I, label="I")
plot!(plt1, t, tot_R, label="R")
plot!(plt1, t, tot_V, label="V")
plot!(plt1, t, tot_N, label="N")
xlabel!(plt1, "time (years)")
ylabel!(plt1, "number")
savefig(plt1, joinpath(figures_dir, "compartments_over_time_$(output_tag).pdf"))

# total mortality over time (after transient)
mort_total = vec(sum(M_mat, dims=2))
sel = t .>= 20.0
plt2 = plot(t[sel], mort_total[sel], label="mortality")
xlabel!(plt2, "time (years)")
ylabel!(plt2, "mortality incidence (per year)")
savefig(plt2, joinpath(figures_dir, "mortality_over_time_$(output_tag).pdf"))

# mean age of infection over time   
mai_over_time = similar(t)
I_div_N_vec   = zeros(n_age)
FOI_vec       = zeros(n_age)

for (k, tk) in enumerate(t)
    S_k  = @view S_mat[k, :]
    I_k  = @view I_mat[k, :]
    R1_k = @view R1_mat[k, :]
    V1_k = @view V1_mat[k, :]
    R2_k = @view R2_mat[k, :]
    V2_k = @view V2_mat[k, :]

    N_k = S_k .+ I_k .+ R1_k .+ R2_k .+ V1_k .+ V2_k .+ 1e-6

    @. I_div_N_vec = I_k / N_k
    mul!(FOI_vec, contact_matrix, I_div_N_vec)

    fsc = 1 + forcing_single * cos(2π * tk)
    @. FOI_vec = β_single * fsc * FOI_vec      # λ(t,a)

    weights = FOI_vec .* S_k                   # λ(t,a) S(t,a)
    denom   = sum(weights) * da
    if denom > 0
        num = sum(age_points .* weights) * da
        mai_over_time[k] = num / denom
    else
        mai_over_time[k] = NaN
    end
end

pltMAI = plot(t, mai_over_time, label="mean age of infection")
xlabel!(pltMAI, "time (years)")
ylabel!(pltMAI, "mean age of infection (years)")
savefig(pltMAI, joinpath(figures_dir, "mean_age_infection_over_time_$(output_tag).pdf"))

# time–age heatmaps
maxI = maximum(I_mat[end, :])
maxM = maximum(M_mat[end, :])

pltS = heatmap(t, age_points, S_mat', xlabel="time (years)", ylabel="age (years)")
savefig(pltS, joinpath(figures_dir, "heatmap_S_$(output_tag).pdf"))

pltI = heatmap(t, age_points, I_mat', xlabel="time (years)", ylabel="age (years)", clim=(0, maxI))
savefig(pltI, joinpath(figures_dir, "heatmap_I_$(output_tag).pdf"))

pltR = heatmap(t, age_points, (R1_mat .+ R2_mat)', xlabel="time (years)", ylabel="age (years)")
savefig(pltR, joinpath(figures_dir, "heatmap_R_$(output_tag).pdf"))

pltV = heatmap(t, age_points, (V1_mat .+ V2_mat)', xlabel="time (years)", ylabel="age (years)")
savefig(pltV, joinpath(figures_dir, "heatmap_V_$(output_tag).pdf"))

pltM = heatmap(t, age_points, M_mat', xlabel="time (years)", ylabel="age (years)", clim=(0, maxM))
savefig(pltM, joinpath(figures_dir, "heatmap_M_$(output_tag).pdf"))

pltN = heatmap(t, age_points, (S_mat .+ I_mat .+ R1_mat .+ R2_mat .+ V1_mat .+ V2_mat)', xlabel="time (years)", ylabel="age (years)")
savefig(pltN, joinpath(figures_dir, "heatmap_N_$(output_tag).pdf"))
 