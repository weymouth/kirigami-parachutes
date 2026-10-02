# Single free-fall run with VTK output of the flow and body motion, for the supplementary videos.
# Not used by any figure. Output: vtk_data/kirigami_N128_H0.25_fall_*.vti (~1.2GB per frame)
include(joinpath(@__DIR__,"..","src","free_fall.jl"))
mkpath(datadir); cd(datadir)

N = 2^7; R = 2N/3.f0

let H = 0.25, θ₀ = 0.4f0
    sim = kirigami_half(N;mem=CuArray,H,fall=true,θ₀); drop!(sim,compute_parameters(N,H,θ₀;R,mem=CuArray))
    X₀(i) = sim.flow.uBC.body.X₀[i]
    motion(a) = (a.flow.f .= 0; a.flow.f[:,:,:,1] .= X₀(1); a.flow.f[:,:,:,2] .= X₀(2); a.flow.f |> Array)
    writer = vtkWriter("kirigami_N$(N)_H$(H)_fall"; attrib=Dict("u"=>vtk_u,"ω"=>vtk_ω,"λ₂"=>vtk_λ₂,"d"=>vtk_d,"motion"=>motion))
    foreach(t->(sim_step!(sim,t); save!(writer,sim)), 0.2:0.05:20.0)
    close(writer)
end
