# Load-dependent mass split: canopy fraction f=1/(1+H/H₁) at fixed total mass, with a point payload at the nose and rotation about the true centre of mass.
# Runs both perturbations for H=1/4,1/2,1,2,4.
# Output: data/load_dependent_check/kirigami_N128_*_fall.jld2, with a summary of the pitch in the log
include(joinpath(@__DIR__,"..","src","free_fall.jl"))
using JLD2
out = joinpath(datadir,"load_dependent_check"); mkpath(out); cd(out)

N = 2^7; R = 2N/3.f0; times = 0.2:0.2:20.0
function run(H,θ₀,dims,name)
    params = compute_parameters(N,H,0.2f0;R,mem=CuArray,T=Float32)
    println("PARAMS $name ",params,"  I/Iₐ=",params.I/params.Iₐ); flush(stdout)
    sim = kirigami_half(N;mem=CuArray,H,fall=true,θ₀,dims); drop!(sim,params)
    data = freefall!(sim,times)
    save_object("kirigami_N$(N)_$(name)_fall.jld2",data)
    println("RESULT $name t_end=",data.t[end],"  θ° every 2R/U: ",round.(rad2deg.(data.θ[10:10:end]),digits=1)); flush(stdout)
end

for H in (0.25f0,0.5f0,1.f0,2.f0,4.f0), θ₀ in (0.2f0,0.4f0)
    run(H,θ₀,(6N,4N,3N÷2),"H$(H)_θ$(θ₀)")
end
