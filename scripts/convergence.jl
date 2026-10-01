# Resolution convergence study of the H=1, 8-ring parachute with prescribed acceleration.
# Output: kirigami_N$(N)_H1_rings8_hist.jld2 and VTK fields, used by figure_3 in figures/figures.jl
include(joinpath(@__DIR__,"..","src","quarter_domain.jl"))
using CUDA,JLD2
mkpath(datadir); cd(datadir)

H = 1; rings = 8; times = 0.05:0.05:3
for N in (64, 96, 128, 192, 256, 384) # N=384 needs ~10GiB of GPU memory
    @show N; flush(stdout)
    sim = kirigami_quarter(N;H,rings,mem=CUDA.CuArray)
    data = drag!(sim,times)
    save_object("kirigami_N$(N)_H$(H)_rings$(rings)_hist.jld2",data)
    writer = vtkWriter("kirigami_N$(N)_H$(H)_rings$(rings)"; attrib=Dict("ω"=>vtk_ω,"λ₂"=>vtk_λ₂,"d"=>vtk_d))
    save!(writer,sim); close(writer)
end
