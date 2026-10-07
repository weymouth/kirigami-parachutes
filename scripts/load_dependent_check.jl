# Check of the load-dependent mass split (canopy fraction f=1-H/H₀) for both perturbations at H=1/4, 1/2 and 1.
# Output: data/load_dependent_check/kirigami_N128_H*_θ*_fall.jld2 and a summary of the pitch peaks in run.log
include(joinpath(@__DIR__,"..","src","free_fall.jl"))
using JLD2
out = joinpath(datadir,"load_dependent_check"); mkpath(out); cd(out)

N = 2^7; R = 2N/3.f0
for H in (0.25f0,0.5f0,1.f0)
    params = compute_parameters(N,H,0.2f0;R,mem=CuArray,T=Float32)
    println("PARAMS H=$H ",params,"  I/Iₐ=",params.I/params.Iₐ); flush(stdout)
    for θ₀ in (0.2f0,0.4f0)
        sim = kirigami_half(N;mem=CuArray,H,fall=true,θ₀,dims=(6N,4N,3N÷2)); drop!(sim,params)
        data = freefall!(sim,0.2:0.2:20.0)
        save_object("kirigami_N128_H$(H)_θ$(θ₀)_fall.jld2",data)
        println("RESULT H=$H θ₀=$θ₀ t_end=",data.t[end],"  θ° every 2R/U: ",round.(rad2deg.(data.θ[10:10:end]),digits=1)); flush(stdout)
    end
end
