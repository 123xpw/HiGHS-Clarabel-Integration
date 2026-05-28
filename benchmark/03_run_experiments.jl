"""
03_run_experiments.jl

Run 24 MPS files × 3 solver configurations:
  HiPO    : solver = simplex       (HiGHS primal simplex, default)
  IPX     : solver = ipm           (HiGHS interior-point / IPX)
  Clarabel: mip_lp_solver = clarabel  (our integration, root-LP only)

Each run: time_limit=3600s, mip_rel_gap=0.01

Results  : benchmark/results.csv
Solver log: benchmark/logs/<case>_<solver>.log

Run:
  julia benchmark/03_run_experiments.jl
  julia benchmark/03_run_experiments.jl --time_limit 7200 --gap 0.005
"""

using Dates

const HIGHS_BIN = abspath(joinpath(@__DIR__, "..", "build_clarabel", "bin", "highs"))
const MPS_DIR   = joinpath(@__DIR__, "mps")
const LOG_DIR   = joinpath(@__DIR__, "logs")
const OPTS_DIR  = joinpath(@__DIR__, "highs_opts")
const RESULTS   = joinpath(@__DIR__, "results.csv")

isfile(HIGHS_BIN) || error(
    "HiGHS binary not found: $HIGHS_BIN\n" *
    "Build: cmake --build build_clarabel -- -j\$(nproc)")

mkpath(LOG_DIR)
mkpath(OPTS_DIR)

# ── CLI args ───────────────────────────────────────────────
time_limit = 3600.0
mip_gap    = 0.01
i = 1
while i <= length(ARGS)
    if ARGS[i] == "--time_limit" && i < length(ARGS)
        time_limit = parse(Float64, ARGS[i+1]); i += 2
    elseif ARGS[i] == "--gap" && i < length(ARGS)
        mip_gap = parse(Float64, ARGS[i+1]); i += 2
    else
        i += 1
    end
end

println("time_limit = $(time_limit)s   mip_rel_gap = $mip_gap")

# ── Write HiGHS options files ──────────────────────────────
# Common options shared by all configs (written via options file so
# mip_lp_solver is also settable — it is not a CLI flag).
common = """
time_limit = $time_limit
mip_rel_gap = $mip_gap
log_to_console = true
"""

configs = [
    (name="IPX",      extra="mip_lp_solver = ipx\n"),
    (name="HiPO",     extra="mip_lp_solver = hipo\n"),
    (name="Clarabel", extra="mip_lp_solver = clarabel\n"),
]

opt_files = Dict{String,String}()
for cfg in configs
    path = joinpath(OPTS_DIR, "$(lowercase(cfg.name)).opt")
    open(path, "w") do f; print(f, common * cfg.extra); end
    opt_files[cfg.name] = path
end

# ── Output parser ──────────────────────────────────────────
function parse_highs_output(text::String)
    obj    = NaN;  time   = NaN
    nodes  = -1;   gap    = NaN
    status = "unknown"

    for line in split(text, '\n')
        # HiGHS 1.14 "Solving report" format
        m = match(r"^\s+Status\s+(.+)", line)
        m !== nothing && (status = strip(m.captures[1]))

        m = match(r"^\s+Primal bound\s+([\d.eE+\-]+)", line)
        m !== nothing && (obj = parse(Float64, m.captures[1]))

        m = match(r"^\s+Timing\s+([\d.eE+\-]+)", line)
        m !== nothing && (time = parse(Float64, m.captures[1]))

        m = match(r"^\s+Nodes\s+(\d+)", line)
        m !== nothing && (nodes = parse(Int, m.captures[1]))

        m = match(r"^\s+Gap\s+([\d.eE+\-]+)%", line)
        m !== nothing && (gap = parse(Float64, m.captures[1]))
    end
    return (obj=obj, time=time, nodes=nodes, gap=gap, status=status)
end

# ── Main loop ──────────────────────────────────────────────
mps_files = sort(filter(f -> endswith(f, ".mps"), readdir(MPS_DIR)))
isempty(mps_files) && error("No MPS files in $MPS_DIR — run 02_generate_mps.jl first")

total_runs = length(mps_files) * length(configs)
run_count  = 0

csv_rows = String["case,month,date,grid,periods,solver,status,objective,time_s,nodes,gap,log_file"]

println("="^72)
println("$(length(mps_files)) MPS × $(length(configs)) solvers = $total_runs runs")
println("="^72)

for mps_file in mps_files
    mps_path = joinpath(MPS_DIR, mps_file)
    base     = replace(mps_file, ".mps" => "")

    # parse grid / date from filename  (e.g. case2383wp_2017-03-14 or case118_2017-01-09)
    m = match(r"^(case\w+)_(\d{4}-\d{2}-\d{2})$", base)
    if m === nothing
        @warn "Unexpected MPS filename: $mps_file"
        continue
    end
    grid, date_str = m.captures
    month   = parse(Int, date_str[6:7])
    periods = startswith(grid, "case2383") ? 72 : 36

    println("\n── $base")

    for cfg in configs
        global run_count += 1
        log_path = joinpath(LOG_DIR, "$(base)_$(lowercase(cfg.name)).log")

        # Skip if already completed (has a feasible result — Optimal or Time limit with valid obj)
        if isfile(log_path) && occursin(r"^\s+Status\s+\S"m, read(log_path, String))
            existing = parse_highs_output(read(log_path, String))
            if existing.status == "Optimal" || (existing.status == "Time limit reached" && !isnan(existing.obj))
                t_str = isnan(existing.time) ? "?" : "$(round(existing.time,digits=2))s"
                println("  [$(lpad(run_count,3))/$(total_runs)] $(rpad(cfg.name, 8))  SKIP (done: $(rpad(existing.status,14)) obj=$(existing.obj)  t=$(t_str))")
                push!(csv_rows,
                    "$base,$month,$date_str,$grid,$periods,$(cfg.name)," *
                    "$(existing.status),$(existing.obj),$(round(coalesce(existing.time,0.0),digits=3))," *
                    "$(existing.nodes),$(existing.gap),$(basename(log_path))")
                continue
            else
                println("  [$(lpad(run_count,3))/$(total_runs)] $(rpad(cfg.name, 8))  RERUN (prev: $(existing.status), no feasible obj)")
            end
        end

        print("  [$(lpad(run_count,3))/$(total_runs)] $(rpad(cfg.name, 8))")
        flush(stdout)

        cmd = `$HIGHS_BIN --model_file $mps_path --options_file $(opt_files[cfg.name])`

        buf    = IOBuffer()
        t_wall = @elapsed begin
            try
                run(pipeline(cmd, stdout=buf, stderr=buf))
            catch
                # HiGHS exits non-zero for INFEASIBLE / TIME_LIMIT — keep captured output
            end
        end
        output = String(take!(buf))

        open(log_path, "w") do f; print(f, output); end

        r = parse_highs_output(output)
        t_str = isnan(r.time) ? "$(round(t_wall,digits=1))s" : "$(round(r.time,digits=2))s"
        println("  $(rpad(r.status, 14)) obj=$(r.obj)  t=$(t_str)  nodes=$(r.nodes)")

        push!(csv_rows,
            "$base,$month,$date_str,$grid,$periods,$(cfg.name)," *
            "$(r.status),$(r.obj),$(round(coalesce(r.time,t_wall),digits=3))," *
            "$(r.nodes),$(r.gap),$(basename(log_path))")
    end
end

# ── Write results ──────────────────────────────────────────
open(RESULTS, "w") do f
    foreach(r -> println(f, r), csv_rows)
end

println()
println("="^72)
println("Results : $RESULTS")
println("Logs    : $LOG_DIR")
println("Runs    : $run_count / $total_runs completed")
println("="^72)
