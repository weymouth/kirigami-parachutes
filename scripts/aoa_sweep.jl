# Prescribed acceleration at a fixed angle of attack θ₀=0.1 for a range of deployment ratios H,
# giving the lift and moment curve slopes.
# Output: kirigami_N128_H$(H)_AoA_fall.jld2, used by figure_4 in figures/figures.jl
include(joinpath(@__DIR__,"..","src","half_domain.jl"))
using JLD2
mkpath(datadir); cd(datadir)

N = 2^7; θ₀ = 0.1f0; R = 2N/3.f0
times = 0.01:0.01:30.0
for H in (0.0,0.25,0.5,1.0,2.0,4.0)
    @show H; flush(stdout)
    sim = kirigami_half(N;mem=CuArray,H,fall=false,θ₀,dims=(6N,3N,3N÷2));
    Xₘ = H==0 ? sim.body.map.x₀ : sim.body.a.b.map.x₀ # moment point in lab frame
    data = drag!(sim,times,R,Xₘ) # run
    save_object("kirigami_N$(N)_H$(H)_AoA_fall.jld2",data)
    writer = vtkWriter("kirigami_N$(N)_H$(H)_AoA_fall"; attrib=Dict("ω"=>vtk_ω,"λ₂"=>vtk_λ₂,"d"=>vtk_d))
    save!(writer,sim); close(writer)
end
