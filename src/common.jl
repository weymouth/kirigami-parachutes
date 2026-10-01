# Code shared by all the simulation and figure scripts
using WaterLily,BiotSavartBCs,StaticArrays,TypedTables
include("paths.jl")

# make a writer with some attributes, need to output to CPU array to save file (|> Array)
using WriteVTK
import WaterLily: @loop,ω,λ₂
vtk_ω(a::AbstractSimulation) = (@loop a.flow.f[I,:] .= ω(I,a.flow.u) over I in inside(a.flow.p); a.flow.f |> Array)
vtk_d(a::AbstractSimulation) = (measure_sdf!(a.flow.σ,a.body,WaterLily.time(a)); a.flow.σ |> Array)
vtk_λ₂(a::AbstractSimulation) = (@inside a.flow.σ[I] = λ₂(I,a.flow.u); a.flow.σ |> Array)
vtk_u(a::AbstractSimulation) = a.flow.u |> Array
