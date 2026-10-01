# Kirigami parachutes

Scripts to reproduce the results of the kirigami parachute paper (`paper/JFM_paper.tex`), simulated with [WaterLily.jl](https://github.com/WaterLily-jl/WaterLily.jl) and the Biot-Savart far-field boundary conditions of [BiotSavartBCs.jl](https://github.com/WaterLily-jl/BiotSavartBCs.jl).

> **Work in progress.** The scripts are being consolidated from the original study; see the open items below.

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

## License

The code in this repository is released under the [MIT license](LICENSE).
