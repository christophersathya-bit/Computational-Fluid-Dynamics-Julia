#!/usr/bin/env bash
# Installs deps, validates all scripts, then runs each one headless. Logs go to ./logs
set -uo pipefail
cd "$(dirname "$0")"

export GKSwstype=100          # headless GR backend (no display needed)
mkdir -p logs
fail=0

command -v julia >/dev/null || { echo "julia not found in PATH"; exit 1; }

echo "== Dependencies"
julia --project=. -e 'using Pkg; "Plots" in keys(Pkg.project().dependencies) || Pkg.add("Plots"); Pkg.instantiate()' \
    || { echo "dependency install failed"; exit 1; }

echo "== Validation"
julia --project=. validate.jl 2>&1 | tee logs/validate.log || fail=1

echo "== Running scripts"
run() {   # run <file> [stdin text]
    local f=$1 inp=${2:-}
    printf -- "-> %s\n" "$f"
    printf "%b" "$inp" | julia --project=. "$f" > "logs/${f%.jl}.log" 2>&1 \
        || { echo "   FAILED (see logs/${f%.jl}.log)"; fail=1; }
}

run Advection_Diffusion_OneDim.jl
run Advection_Diffusion_OneDim_Upwind.jl
run Demo_1DProblem.jl
run TwoDim_Conv_Diff.jl
run DemoScript.jl          "1.5\n0.1\n"
run Programming_Demo.jl    "1.5\n0.1\n"
run DemoScript_2.jl
run Programming_Demo_2.jl  "10\n"

if [ "$fail" -eq 0 ]; then echo "== ALL OK"; else echo "== SOME STEPS FAILED (see logs/)"; fi
exit "$fail"
