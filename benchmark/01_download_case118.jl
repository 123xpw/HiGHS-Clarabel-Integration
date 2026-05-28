"""
01_download_case118.jl

Download 1 random day per month (2017) for case118 into
  testdata/UnitCommitment_Data/case118/

January  : copy a random day from instances/matpower/case118/ (31 days exist)
Feb–Dec  : download from axavier.org/UnitCommitment.jl/0.3

Run from anywhere:
  julia .../highs-clarabel/benchmark/01_download_case118.jl
"""

import Downloads
using Dates
using Random

const SEED        = 42
const BASE_URL    = "https://axavier.org/UnitCommitment.jl/0.3/instances/matpower"
const UC_ROOT     = abspath(joinpath(@__DIR__, "..", "..", "..", "UnitCommitment.jl"))
const MATPOWER118 = joinpath(UC_ROOT, "instances", "matpower", "case118")
const OUT_DIR     = joinpath(UC_ROOT, "testdata", "UnitCommitment_Data", "case118")

mkpath(OUT_DIR)

rng = MersenneTwister(SEED)

println("="^60)
println("Download / copy case118 – one day per month (2017)")
println("UC root : $UC_ROOT")
println("Output  : $OUT_DIR")
println("="^60)

selected = Dict{Int,String}()   # month => "YYYY-MM-DD"

# ── January: copy from local matpower data ─────────────────
jan_files = sort(filter(f -> startswith(f, "2017-01-") && endswith(f, ".json.gz"),
                        readdir(MATPOWER118)))
isempty(jan_files) && error("No January files found in $MATPOWER118")
jan_pick  = rand(rng, jan_files)
src = joinpath(MATPOWER118, jan_pick)
dst = joinpath(OUT_DIR, jan_pick)
if !isfile(dst)
    cp(src, dst)
    println("Jan  COPIED : $jan_pick")
else
    println("Jan  exists : $jan_pick")
end
selected[1] = replace(jan_pick, ".json.gz" => "")

# ── Feb–Dec: download from axavier.org ─────────────────────
for month in 2:12
    days_in_month = Dates.daysinmonth(Date(2017, month, 1))
    day       = rand(rng, 1:(days_in_month - 1))   # avoid last day → D+1 always exists
    date_str  = Dates.format(Date(2017, month, day), "yyyy-mm-dd")
    filename  = "$date_str.json.gz"
    filepath  = joinpath(OUT_DIR, filename)
    tag       = lpad(Dates.monthabbr(month), 3)

    if isfile(filepath)
        println("$tag  exists : $filename")
        selected[month] = date_str
        continue
    end

    url = "$BASE_URL/case118/$filename"
    print("$tag  downloading $filename ... ")
    flush(stdout)
    try
        Downloads.download(url, filepath)
        println("OK")
        selected[month] = date_str
    catch e
        println("FAILED – $(sprint(showerror, e))")
        println("    The server may not have data for this date.")
        println("    Consider running generate_year_cases.jl as fallback.")
    end
end

println()
println("Selected dates:")
for m in 1:12
    tag = haskey(selected, m) ? selected[m] : "MISSING"
    println("  $(lpad(m,2)) $(rpad(Dates.monthname(m),9)): $tag")
end
println()
println("Done. Files in: $OUT_DIR")
