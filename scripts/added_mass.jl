# Linear and rotational added mass of the 16-ring parachute for a range of deployment ratios H.
# Output: kirigami_parameters.jld2, used by figure_4 in figures/figures.jl and by the free-fall runs
include(joinpath(@__DIR__,"..","src","free_fall.jl"))
using JLD2
mkpath(datadir); cd(datadir)

N = 2^7; R = 2N/3.f0
params_all = []
for H in (0.0,0.25,0.5,1.0,2.0,4.0)
    @show H; flush(stdout)
    params = compute_parameters(N,H,0;R,mem=CuArray,T=Float32)
    push!(params_all, params)
end; save_object("kirigami_parameters.jld2",params_all)
