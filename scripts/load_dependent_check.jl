# Load-dependent mass split: canopy fraction f=1/(1+H/H₁) at fixed total mass, with a point payload at the nose and rotation about the true centre of mass.
# Runs the two larger domains of the H=1, θ₀=0.2 domain sweep (the H and θ₀ sweep is in git history).
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

for dims in ((6N,4N,3N),(8N,4N,3N))
    run(1.f0,0.2f0,dims,"$(dims[1])x$(dims[2])x$(dims[3])")
end
