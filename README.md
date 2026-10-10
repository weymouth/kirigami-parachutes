# Kirigami parachutes

Scripts for the kirigami parachute paper (`paper/JFM_paper.tex`), using [WaterLily.jl](https://github.com/WaterLily-jl/WaterLily.jl), [BiotSavartBCs.jl](https://github.com/WaterLily-jl/BiotSavartBCs.jl) and [WaterLilyNarrowBand.jl](https://github.com/weymouth/WaterLilyNarrowBand.jl).

> **Still being finalized.** WaterLily and WaterLilyNarrowBand come from GitHub until their next releases (`[sources]` in `Project.toml`).

## Setup

Requires Julia `v1.11` and an NVIDIA GPU.
```bash
git clone https://github.com/weymouth/kirigami-parachutes
cd kirigami-parachutes
julia --project -e "using Pkg; Pkg.instantiate()"
```
Without a display, precompile GLMakie with `xvfb-run`. `Manifest.toml` is committed once the paper is accepted.

## Reproducing the paper

Simulations write to `data/`, figures to `paper/fig`:
```bash
julia --project -t auto scripts/convergence.jl
julia --project -t auto figures/figures.jl
```

| Simulation script | Data | Used in |
|---|---|---|
| `scripts/convergence.jl` | `kirigami_N*_H1_rings8_hist.jld2`, VTK fields | `figure_3` |
| `scripts/rings_sweep.jl` | `kirigami_N256_H*_rings*_hist.jld2` | `figure_4` |
| `scripts/deployment_sweep.jl` | `kirigami_N256_H*_hist.jld2` | `figure_4` |
| `scripts/added_mass.jl` | `kirigami_parameters.jld2` | `figure_4`, free fall |
| `scripts/aoa_sweep.jl` | `kirigami_N128_H*_AoA_fall.jld2` | `figure_4` |
| `scripts/free_fall.jl` | `kirigami_N128_*_fall.jld2` | `figure_5`, `figure_6` |
| `scripts/load_dependent_check.jl` | `load_dependent_check/kirigami_N128_*_fall.jld2` | load-dependent mass check |
| `scripts/free_fall_video.jl` | VTK fields (~500 GB) | supplementary videos |
| `scripts/sphere.jl` | `sphere_*.jld2` | `figure_A2` |
| `scripts/impulsive_circle.jl` | `ImpCircle_results.jld2` | `figure_A1` |
| `scripts/benchmark.jl`, `scripts/narrowband_timing.jl` | log only | performance checks |

Fig. 1 is a ParaView rendering of the `deployment_sweep.jl` VTK output and Fig. 2 an Inkscape schematic. `figure_A2` uses GLMakie and needs a display.

Axial runs use a quarter domain (`src/quarter_domain.jl`), runs with motion in x-y a half domain (`src/half_domain.jl`); both build the parachute in `src/kirigami.jl`. Free-fall dynamics are in `src/free_fall.jl`.

## License

[MIT](LICENSE)
