# 7_1_CertifiedSignSequence.jl — README

## What this script does

`7_1_CertifiedSignSequence.jl` generates **rigorously certified sign sequences** (`L`/`R`) for the Lorenz system by detecting **upward crossings** of the plane `z = 27`. Each certified crossing is labeled by the **sign of `y`** on a rigorously validated time bracket around the crossing (`L` if `y<0`, `R` if `y>0`).

The run has two phases:

1. **Fast, uncertified burn-in** (e.g., RK4 via OrdinaryDiffEq) from a random initial point to land near the strange attractor.
2. **Certified, chunked integration** using Taylor models to:
   - propagate a validated flowpipe,
   - over-approximate each segment to a box,
   - **bracket** each candidate crossing with local re-integration and **bisection** until the time bracket is ≤ a target width,
   - **prove** directionality via an interval bound `inf(dz/dt) > 0` on the hull of the start/mid/end boxes,
   - classify the symbol by the uniform sign of `y` on the bracket.

The script **stops at the first crossing it cannot certify**, by design, so you can see where precision becomes insufficient under your chosen settings.

It then **optionally starts another sequence** by taking a point estimate (box midpoint) at the stopping location, doing another short **uncertified hop**, and resuming certified integration—repeating until the requested number of sequences is produced. A plain-text report is written to disk (seed, parameters, initial/burn-in points, and the sign sequences).

---

## Key packages (for certification)

- **ReachabilityAnalysis.jl** — validated reachability for ODEs; the script uses the **TMJets** Taylor-model algorithm for rigorous flowpipes and local re-integration.
- **LazySets.jl** — set representations; we **over-approximate** flowpipe segments to `Hyperrectangle`s to get per-coordinate bounds.
- **IntervalArithmetic.jl** — outward-rounded intervals; used to bound `x`, `y`, `z`, and the **derivative** `dz/dt = x*y − (8/3)z` on brackets for directionality.
- **OrdinaryDiffEq.jl** — fast (non-rigorous) integrators for **burn-in** and random hops to re-seed near the attractor.

> Notes on rigor: Certification relies on (i) a sign change of `z − 27` between validated boxes with **strict** `sup(z_lo) < 27 < inf(z_hi)`, (ii) **positive** interval lower bound on `dz/dt` across the **hull** of endpoint/midpoint boxes, and (iii) **uniform** sign of `y` on the accepted time bracket.

---

## Requirements

- Julia **1.11** (or newer in the 1.x series)
- Packages:
  ```julia
  import Pkg
  Pkg.activate()  # global env
  Pkg.add([
      "ReachabilityAnalysis",
      "LazySets",
      "IntervalArithmetic",
      "OrdinaryDiffEq",
  ])

*Ignore TaylorIntegration package errors*

## Running the script
## *You will get error messages!* 
Some (all?) are related to the TaylorIntegration package. Ignore these. If you see something like this below the error message, everything is working.

```terminal
Starting run at 2025-09-07T13:46:26.557. Version v7.3
Requested sequences: 2  |  burn-in: 5.0 s  |  output: CertifiedSignSequence_report.txt
Random initial point: (-13.975531465029519, -8.139724144449527, 24.980481884986823)
Post burn-in point 1 = (-10.20507814064993, -11.042315229746569, 28.12924724332252)

===== Generating certified sign sequence 1 =====
```

How to run:

```bash
julia /path/to/7_1_CertifiedSignSequence.jl [flags...]
```

Common flags

--num-sequences=<Int> : how many sign sequences to produce.

--seed=<Int> : single seed for reproducible randomness (initial point and burn-in hops).

--burn-in=<Float> : burn-in time in seconds (uncertified).

--output=<String> : output filename for the text report.

--target-time-width=<Float> : target bracket width (seconds) in the local bisection (smaller = stricter, slower).

--global-t-max=<Float> : overall time budget for certified integration.

--chunk-len=<Float> : certified integration chunk size (seconds) before scanning and continuing.

--hop-min=<Float> / --hop-max=<Float> : min/max seconds for uncertified random hops between sequences.

--eps-box=<Float> : tiny box radius around point estimates when switching from uncertified to certified phases.

--init-range=a,b,c,d,e,f : initial random box ranges (xl,xh,yl,yh,zl,zh).

All flags are optional; sane defaults are provided. Use --seed for reproducible runs (the script seeds all random draws from it).

Examples:

Produce two sign sequences:

```bash
julia "/path/to/7_1CertifiedSignSequence.jl" --num-sequences=2
```

Fixed seed, shorter brackets, custom output file

```bash
julia "/path/to/7_1_CertifiedSignSequence.jl" \
  --num-sequences=2 \
  --seed=930415324 \
  --target-time-width=0.05 \
  --output="CertifiedSignSequences_run.txt"
```

Output

*Ignore

A .txt report is written with:

Script version
Program run time
Random seed (ensures full reproducibility of initial point and hops)
All parameter settings
Initial random point and post burn-in points
Sign sequences, one per line at the end