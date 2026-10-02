# Kirigami parachutes

Scripts to reproduce the results of the kirigami parachute paper (`paper/JFM_paper.tex`), simulated with [WaterLily.jl](https://github.com/WaterLily-jl/WaterLily.jl) and the Biot-Savart far-field boundary conditions of [BiotSavartBCs.jl](https://github.com/WaterLily-jl/BiotSavartBCs.jl).

> **Still being finalized.** This repository is being consolidated from the original study and checked against the published results, so scripts and data formats may still change. The environment uses WaterLily 1.9.1 and, until its registration completes, pins BiotSavartBCs to a commit on GitHub.

## Layout

| Folder | Contents |
|---|---|
| `paper/` | LaTeX source of the paper, reviewer responses, and figures (`paper/fig`) |
| `scripts/` | Simulation scripts that generate the data |
| `figures/` | Scripts that turn the data into the figures in `paper/fig` |
| `src/` | Code shared by the scripts |
| `data/` | Simulation output (not version controlled) |

The simulation data is not stored in this repository; running the scripts regenerates it.

## Setup

Requires Julia `v1.11` or later and an NVIDIA GPU for the 3D simulations.
```bash
git clone https://github.com/weymouth/kirigami-parachutes
cd kirigami-parachutes
julia --project -e "using Pkg; Pkg.instantiate()"
```
The committed `Manifest.toml` pins the exact package versions used for the paper.

## Reproducing the paper

Run the simulation scripts, which write their output to `data/`, e.g.
```bash
julia --project -t auto scripts/convergence.jl
```
and then the figure scripts, which write to `paper/fig`:
```bash
julia --project figures/figures.jl
```

| Simulation script | Data | Used in |
|---|---|---|
| `scripts/convergence.jl` | `kirigami_N*_H1_rings8_hist.jld2`, VTK fields | `kirigami_convergence.png` (`figure_3`) |
| `scripts/rings_sweep.jl` | `kirigami_N256_H*_rings*_hist.jld2` | `kirigami_Cd_time.png` (`figure_4`) |
| `scripts/deployment_sweep.jl` | `kirigami_N256_H*_hist.jld2` | `kirigami_Cd_time.png` (`figure_4`) |
| `scripts/added_mass.jl` | `kirigami_parameters.jld2` | `kirigami_Cd_time.png` (`figure_4`), free-fall runs |
| `scripts/aoa_sweep.jl` | `kirigami_N128_H*_AoA_fall.jld2` | `kirigami_Cd_time.png` (`figure_4`) |
| `scripts/free_fall.jl` | `kirigami_N128_*_fall.jld2` | `kirigami_domain_sweep.png` (`figure_5`), `kirigami_results.png` (`figure_6`) |
| `scripts/sphere.jl` | `sphere_*.jld2`, `sphere_320x128x128_t53.jld2` (flow state) | `validation_sphere.png` (`figure_A2`) |
| `scripts/impulsive_circle.jl` | `ImpCircle_results.jld2` | `ImpCircle_results.png` (`figure_A1`) |

`figure_A2` also renders the sphere wake from the saved flow state (`sphere3_zoom.png`, embedded in `validation_sphere.png`) with WaterLily's Makie `viz!`, which needs a display for GLMakie. The remaining figures are hand-made: `flow_render_2.png` (Fig. 1) is a ParaView rendering of the λ₂ isosurfaces in the VTK output of `scripts/deployment_sweep.jl`, and `multilevel_domain.svg` (Fig. 2) is a schematic drawn in Inkscape.

Simulations with purely axial motion (convergence, rings and deployment sweeps) use a quarter domain with symmetry planes in y and z (`src/quarter_domain.jl`). Simulations with motion in the x-y plane (added mass, angle of attack and free fall) use a half domain with a symmetry plane in z (`src/half_domain.jl`). Each overrides the Biot-Savart symmetry function, so only include one of the two in a Julia session. Both build the parachute from `src/kirigami.jl`. The free-fall added-mass calculation and dynamics are in `src/free_fall.jl`: the falling parachute has mass m=1.5ρR³ concentrated at its nose, so its own moment of inertia about the nose is zero and the rotational inertia is the measured added inertia.

## License

The code in this repository is released under the [MIT license](LICENSE).
