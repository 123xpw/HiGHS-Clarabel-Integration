"""
02_generate_mps.jl

Generate 24 MPS files for the Clarabel benchmark:
  - 12 × case2383wp  (72-period, power-trajectory formulation)
  - 12 × case118     (36-period, standard formulation)

One random day per month is selected from 2017 data (seed=42).
  case2383wp : testdata/UnitCommitment_Data/case2383wp/  (365 days)
  case118    : testdata/UnitCommitment_Data/case118/     (12 days, after 01_download)

Output MPS files : highs-clarabel/benchmark/mps/
Selected dates   : highs-clarabel/benchmark/selected_dates.csv

Run from UnitCommitment.jl root:
  julia --project=. .../highs-clarabel/benchmark/02_generate_mps.jl
Or with explicit project:
  julia .../highs-clarabel/benchmark/02_generate_mps.jl
"""

using Pkg
const UC_ROOT = abspath(joinpath(@__DIR__, "..", "..", "..", "UnitCommitment.jl"))

# Use this benchmark's own Project.toml (has JuMP declared).
# Pkg.develop links the local UnitCommitment source.
Pkg.activate(@__DIR__)
Pkg.develop(path=UC_ROOT)
Pkg.instantiate()

using UnitCommitment
using JuMP
using Dates
using Random

const SEED          = 42
const MPS_DIR       = joinpath(@__DIR__, "mps")
const BENCH_DIR     = @__DIR__

# power-trajectory ramp parameters (expressed in 40-min base time steps)
const BASE_UD_36    = 2
const BASE_DD_36    = 2
const BASE_TIME_MIN = 40

mkpath(MPS_DIR)

# ── Helpers ────────────────────────────────────────────────

"""Pick 1 file per calendar month from a directory of YYYY-MM-DD.json.gz files."""
function pick_one_per_month(data_dir::String, rng::AbstractRNG; avoid_last_day::Bool=false)
    files = sort(filter(f -> occursin(r"^\d{4}-\d{2}-\d{2}\.json\.gz$", f), readdir(data_dir)))
    by_month = Dict{Int, Vector{String}}()
    for f in files
        m = parse(Int, f[6:7])
        push!(get!(by_month, m, String[]), f)
    end
    selected = Dict{Int, String}()
    for month in 1:12
        pool = get(by_month, month, String[])
        isempty(pool) && continue
        if avoid_last_day
            # keep only days where day < days_in_month so D+1 exists in the dataset
            pool = filter(pool) do f
                d = parse(Int, f[9:10])
                d < Dates.daysinmonth(Date(2017, month, 1))
            end
        end
        isempty(pool) && continue
        selected[month] = rand(rng, pool)
    end
    return selected
end

"""Build startup/shutdown trajectory curves from generator limits."""
function build_curves(instance; base_ud=BASE_UD_36, base_dd=BASE_DD_36, base_dt=BASE_TIME_MIN)
    sc = instance.scenarios[1]
    dt = sc.time_step
    su_steps = div(base_ud * base_dt, dt)
    sd_steps = div(base_dd * base_dt, dt)
    (su_steps > 0 && sd_steps > 0) || error("Steps must be positive (dt=$dt)")

    curves = Dict{String, NamedTuple{(:startup,:shutdown),Tuple{Vector{Float64},Vector{Float64}}}}()
    for g in sc.thermal_units
        (g.startup_limit > 0.0 && g.shutdown_limit > 0.0) || continue
        curves[g.name] = (
            startup  = [g.startup_limit  * i / su_steps for i in 1:su_steps],
            shutdown = [g.shutdown_limit * (sd_steps-i+1) / sd_steps for i in 1:sd_steps],
        )
    end
    @info "Curves: $(length(curves)) units | dt=$(dt)min | su=$(su_steps)steps | sd=$(sd_steps)steps"
    return curves
end

# ── Main ───────────────────────────────────────────────────

rng = MersenneTwister(SEED)
csv_rows = ["case,month,date,periods,pt,mps_file"]

println("="^70)
println("Generating MPS files  (seed=$SEED)")
println("UC root : $UC_ROOT")
println("MPS out : $MPS_DIR")
println("="^70)

# ─── 1. case2383wp : 72-period + power-trajectory ──────────────────

println("\n[1/2] case2383wp  (72-period + power-trajectory)")
wp_dir = joinpath(UC_ROOT, "testdata", "UnitCommitment_Data", "case2383wp")
isdir(wp_dir) || error("Not found: $wp_dir")

wp_sel = pick_one_per_month(wp_dir, rng; avoid_last_day=true)
missing_m = setdiff(1:12, keys(wp_sel))
isempty(missing_m) || @warn "case2383wp missing months: $missing_m"

for month in sort(collect(keys(wp_sel)))
    fname    = wp_sel[month]
    date_str = replace(fname, ".json.gz" => "")
    cur_path = joinpath(wp_dir, fname)

    next_date = Dates.format(Date(date_str) + Day(1), "yyyy-mm-dd")
    nxt_path  = joinpath(wp_dir, "$next_date.json.gz")
    isfile(nxt_path) || error("Next-day file missing: $nxt_path")

    print("  month $(lpad(month,2)) | $date_str | 72T+PT ... ")
    flush(stdout)

    inst_cur = UnitCommitment.read(cur_path)
    inst_nxt = UnitCommitment.read(nxt_path)
    UnitCommitment.clear_power_trajectories!(inst_cur)

    inst72 = UnitCommitment.convert_to_subhourly(inst_cur, inst_nxt)
    curves = build_curves(inst72)

    has_pt = !isempty(curves)
    if has_pt
        UnitCommitment.set_power_trajectories!(inst72, curves)
        formulation = UnitCommitment.Formulation(
            power_trajectories = UnitCommitment.xxx2005.PowerTrajectories(),
        )
    else
        @warn "No PT curves for $date_str (all startup_limit=0?); using standard formulation"
        formulation = UnitCommitment.Formulation()
    end

    mps_base = "case2383wp_$date_str"
    mps_path = joinpath(MPS_DIR, "$mps_base.mps")

    model = UnitCommitment.build_model(
        instance       = inst72,
        formulation    = formulation,
        variable_names = true,
    )
    JuMP.write_to_file(model, mps_path)
    sz = round(filesize(mps_path)/1024/1024, digits=1)
    println("OK  ($sz MB)")

    push!(csv_rows, "case2383wp,$month,$date_str,72,$has_pt,$mps_base.mps")
end

# ─── 2. case118 : 36-period standard ───────────────────────────────

println("\n[2/2] case118  (36-period, standard)")
c118_dir = joinpath(UC_ROOT, "testdata", "UnitCommitment_Data", "case118")
isdir(c118_dir) || error("Not found: $c118_dir — run 01_download_case118.jl first")

c118_sel = pick_one_per_month(c118_dir, rng)
missing_m = setdiff(1:12, keys(c118_sel))
isempty(missing_m) || @warn "case118 missing months: $missing_m"

for month in sort(collect(keys(c118_sel)))
    fname    = c118_sel[month]
    date_str = replace(fname, ".json.gz" => "")
    cur_path = joinpath(c118_dir, fname)

    print("  month $(lpad(month,2)) | $date_str | 36T    ... ")
    flush(stdout)

    inst = UnitCommitment.read(cur_path)

    mps_base = "case118_$date_str"
    mps_path = joinpath(MPS_DIR, "$mps_base.mps")

    model = UnitCommitment.build_model(instance=inst, variable_names=true)
    JuMP.write_to_file(model, mps_path)
    sz = round(filesize(mps_path)/1024/1024, digits=2)
    println("OK  ($sz MB)")

    push!(csv_rows, "case118,$month,$date_str,36,false,$mps_base.mps")
end

# ─── Save selected dates ────────────────────────────────────────────

csv_path = joinpath(BENCH_DIR, "selected_dates.csv")
open(csv_path, "w") do f
    foreach(r -> println(f, r), csv_rows)
end

println()
println("="^70)
println("MPS files : $MPS_DIR  ($(length(csv_rows)-1) files)")
println("Date log  : $csv_path")
println("="^70)
