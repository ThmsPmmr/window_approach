# Window matheuristic for the single row facility layout problem (SRFLP)

Code and benchmark instances for:
> T. Pammer, M. Sinnl (2026). A matheuristic for solving the single row facility layout problem. *Computers & Operations Research*. <https://doi.org/10.1016/j.cor.2026.107397>

## Requirements
- Julia 1.11 (all package versions pinned in `Manifest.toml`)
- Gurobi 12.0.0 with a valid licence (free for academics)

## Installation
    julia --project=. -e "using Pkg; Pkg.instantiate()"

## Run one instance
    julia --project=. src/window_matheuristic.jl <instance> <cooling_factor> <p_swap> <inner_runs_factor> <acceleration 0|1> [time_limit] [seed]
    # example with the paper settings:
    julia --project=. src/window_matheuristic.jl instances/Anjos-large/SRFLP_Rubio_AnKeVa_200_01.txt 0.976 0.323 87 1
If `time_limit` is omitted or `0`, the time limits from the paper are used: 3600 s (n ≤ 300), 7200 s (n ≤ 500), 10800 s (n ≤ 1000).
For a quick test, use a short limit (it must be larger than n/2 seconds, otherwise the simulated annealing phase is skipped), e.g. `... 0.976 0.323 87 1 240`.
The last line of the output is a CSV-style summary (instance, parameters, objective value, ...).

## Running all experiments
`scripts/run_parallel.sh` (Linux/macOS, requires [GNU parallel](https://www.gnu.org/software/parallel/)) runs all instances
of a folder in parallel, one run per instance per call:

    scripts/run_parallel.sh instances/Anjos-large anjos 0.976 0.323 87 1 0 32

`0` as the time limit uses the limits from the paper, and `32` is the number of parallel jobs.
Each call appends one run to `output/anjos/<instance>.log`; for the paper it was called 10 times per instance set.
Running all experiments takes about 1,300 CPU hours in total.

## Experimental setup
Each instance was solved in 10 independent runs with cooling factor 0.976, swap probability 0.323,
inner runs factor 87 (× n) and acceleration on, each run on a single core of an AMD EPYC 7543 with 6 GB of memory.
The parameter values were tuned with irace (see Section 5.2 of the paper).
Because the runs are time-limited and randomised, re-runs will differ slightly from the reported values.

## Instances
Anjos-large, Sko-large and Palubeckis-large come from <original papers>. Please cite those sources if you use them.

## License
Code: MIT. The benchmark instances belong to their original authors.

