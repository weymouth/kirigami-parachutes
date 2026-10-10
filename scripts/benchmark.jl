# Time the kirigami measurement step on WaterLily master (b6cf648) and the fused force PR (8a3322e), with and without NarrowBand.
# Usage: julia scripts/benchmark.jl [N]. Makes one environment per WaterLily version in data/.
if length(ARGS) < 2
    include(joinpath(@__DIR__,"..","src","paths.jl"))
    using Pkg
    N = isempty(ARGS) ? "128" : ARGS[1]
    for (name,rev) in ("master"=>"b6cf648414f5808e9ddc26675505bb07523e9ee9", "fused"=>"8a3322ef7f00e2df9db06ef94a612125baba3eb1")
        env = joinpath(datadir,"benchmark_$(name)_$(rev[1:7])")
        if !isfile(joinpath(env,"Manifest.toml"))
            Pkg.activate(env)
            Pkg.add([PackageSpec(url="https://github.com/WaterLily-jl/WaterLily.jl",rev=rev),
                     PackageSpec(name="BiotSavartBCs",version="1"),
                     PackageSpec(url="https://github.com/weymouth/WaterLilyNarrowBand.jl"),
                     "CUDA","StaticArrays","TypedTables","WriteVTK"])
        end
        run(`$(Base.julia_cmd()) -t auto --project=$env $(@__FILE__) $N $name`)
    end
    exit()
end

include(joinpath(@__DIR__,"..","src","half_domain.jl"))
using WaterLilyNarrowBand, Printf
N, name = parse(Int,ARGS[1]), ARGS[2]

body_map(b::NarrowBand) = body_map(b.body)
force_moment(x,sim) = isdefined(WaterLily,:total_force_and_moment) ? WaterLily.total_force_and_moment(x,sim) :
                      (WaterLily.total_force(sim), WaterLily.pressure_moment(x,sim))
timed(f) = (CUDA.synchronize(); @elapsed (f(); CUDA.synchronize()))

for narrow in (false,true)
    sim = kirigami_half(N;mem=CuArray,H=1,θ₀=0.2f0,dims=(6N,4N,3N÷2))
    narrow && (sim.body = NarrowBand(sim.body,size(sim.flow.p).-2;interior=false,mem=CuArray))
    θ, t = 0.2f0, zeros(3)
    for n in 1:25
        x = pivot(body_map(sim.body))
        f = timed(()->force_moment(x,sim))
        θ += 0.005f0; sim.body = setmap(sim.body;θ=SA[0,0,θ])
        m = timed(()->measure!(sim))
        s = timed(()->sim_step!(sim;remeasure=false))
        n > 5 && (t .+= (f,m,s)./20)
    end
    @printf("%-6s %-10s forces %.4fs  measure! %.4fs  flow %.4fs  total %.4fs\n",
            name,narrow ? "NarrowBand" : "SetBody",t...,sum(t))
    sim = nothing; GC.gc(true); CUDA.reclaim()
end
