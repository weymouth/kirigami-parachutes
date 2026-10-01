# Deployment ratio sweep of the 16-ring parachute with prescribed acceleration,
# including negative H to check the symmetry of the results.
# Output: kirigami_N256_H$(H)_hist.jld2, used by figure_4 in figures/figures.jl
include(joinpath(@__DIR__,"..","src","quarter_domain.jl"))
using CUDA,JLD2
mkpath(datadir); cd(datadir)

N = 2^8; times = 0.05:0.05:3
Hs = 0.5 .^ (-2:2)
Hs = [-Hs; 0; reverse(Hs)] # include negative H for checking symmetry
for H ∈ Hs
    @show H; flush(stdout)
    sim = kirigami_quarter(N;H,mem=CUDA.CuArray)
    data = drag!(sim,times)
    save_object("kirigami_N$(N)_H$(H)_hist.jld2",data)
    writer = vtkWriter("kirigami_N$(N)_H$(H)"; attrib=Dict("ω"=>vtk_ω,"λ₂"=>vtk_λ₂,"d"=>vtk_d))
    save!(writer,sim); close(writer)
end
