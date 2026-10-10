# Time the kirigami measurement step on WaterLily master (b6cf648) and setbody-sdf (on the fused force PR), with and without NarrowBand.
# Usage: julia scripts/benchmark.jl [N]. Makes one environment per WaterLily version in data/.
if length(ARGS) < 2
    include(joinpath(@__DIR__,"..","src","paths.jl"))
    using Pkg
    N = isempty(ARGS) ? "128" : ARGS[1]
    for (name,rev) in ("master"=>"b6cf648414f5808e9ddc26675505bb07523e9ee9", "sdf"=>"setbody-sdf")
        env = joinpath(datadir,"benchmark_$(name)_$(rev[1:7])")
        if !isfile(joinpath(env,"Manifest.toml"))
            Pkg.activate(env)
            Pkg.add([PackageSpec(url="https://github.com/WaterLily-jl/WaterLily.jl",rev=rev),
                     PackageSpec(name="BiotSavartBCs",version="1"),
                     PackageSpec(url="https://github.com/weymouth/WaterLilyNarrowBand.jl"),
                     PackageSpec.(["CUDA","StaticArrays","TypedTables","WriteVTK"])...])
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
println("$name on ",CUDA.name(CUDA.device()))

for narrow in (false,true)
    sim = kirigami_half(N;mem=CuArray,H=1,θ₀=0.2f0,dims=(6N,4N,3N÷2))
    narrow && (sim.body = NarrowBand(sim.body,size(sim.flow.p).-2;interior=false,mem=CuArray))
    θ, t = 0.2f0, zeros(5)
    for n in 1:25
        x = pivot(body_map(sim.body))
        f = timed(()->force_moment(x,sim))
        θ += 0.005f0; sim.body = setmap(sim.body;θ=SA[0,0,θ])
        tᵢ = sum(sim.flow.Δt)
        m = timed(()->measure!(sim.flow,sim.body;t=tᵢ,ϵ=sim.ϵ))
        d = timed(()->measure_sdf!(sim.flow.σ,sim.body,tᵢ;fastd²=(2+sim.ϵ)^2)) # again, to time it alone
        p = timed(()->WaterLily.update!(sim.pois))
        s = timed(()->sim_step!(sim;remeasure=false))
        n > 5 && (t .+= (f,d,m,p,s)./20)
    end
    f,d,m,p,s = t
    @printf("%-6s %-10s forces %.4fs  measure_sdf! %.4fs  measure!(flow,body) %.4fs  update!(pois) %.4fs  flow %.4fs  total %.4fs",
            name,narrow ? "NarrowBand" : "SetBody",f,d,m,p,s,f+m+p+s)
    narrow ? @printf("  active cells %d of %d\n",count(sim.body.active),length(sim.body.active)) : println()
    sim = nothing; GC.gc(true); CUDA.reclaim()
end
