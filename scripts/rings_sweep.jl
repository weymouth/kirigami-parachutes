# Number-of-rings sweep with prescribed acceleration, for several deployment ratios H.
# Output: kirigami_N256_H$(H)_rings$(rings)_hist.jld2, used by figure_4 in figures/figures.jl
include(joinpath(@__DIR__,"..","src","quarter_domain.jl"))
using CUDA,JLD2
mkpath(datadir); cd(datadir)

N = 2^8; times = 0.05:0.05:3
for H in (0.5,1,2,4), rings in 4:4:20
    @show H,rings; flush(stdout)
    sim = kirigami_quarter(N;H,rings,mem=CUDA.CuArray)
    data = drag!(sim,times)
    save_object("kirigami_N$(N)_H$(H)_rings$(rings)_hist.jld2",data)
    writer = vtkWriter("kirigami_N$(N)_H$(H)_rings$(rings)"; attrib=Dict("ω"=>vtk_ω,"λ₂"=>vtk_λ₂,"d"=>vtk_d))
    save!(writer,sim); close(writer)
end
