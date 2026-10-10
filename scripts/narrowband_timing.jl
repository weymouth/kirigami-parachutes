# Time the pieces of NarrowBand's measure_sdf! on the kirigami. Run in a benchmark environment:
# julia --project=data/benchmark_sdf_setbody scripts/narrowband_timing.jl [N]
include(joinpath(@__DIR__,"..","src","half_domain.jl"))
using WaterLilyNarrowBand, Printf
import WaterLilyNarrowBand: active_sdf, touching!, activate!
N = isempty(ARGS) ? 128 : parse(Int,ARGS[1])
timed(f) = (f(); CUDA.synchronize(); minimum(@elapsed((f(); CUDA.synchronize())) for _ in 1:10))

sim = kirigami_half(N;mem=CuArray,H=1,θ₀=0.2f0,dims=(6N,4N,3N÷2))
a = NarrowBand(sim.body,size(sim.flow.p).-2;interior=false,mem=CuArray)
d, fd², t = sim.flow.σ, (2+sim.ϵ)^2, 0f0
for θ in 0.2f0 .+ 0.005f0*(1:5)
    global a = setmap(a;θ=SA[0,0,θ]); measure_sdf!(d,a,t;fastd²=fd²)
end
body, active, farinside, near = a.body, a.active, a.farinside, a.near

results = ["measure_sdf!"           => timed(()->measure_sdf!(d,a,t;fastd²=fd²)),
           "  all(active)"          => timed(()->all(active)),
           "  narrow kernel"        => timed(()->(@inside d[I] = active_sdf(body,active,I,t,fd²))),
           "  touching!"            => timed(()->touching!(farinside,d,active,fd²)),
           "  exact cells + any"    => timed(()->(@inside farinside[I] = d[I]^2<fd²; any(farinside))),
           "  activate!"            => timed(()->activate!(a,d,fd²)),
           "one Bool @inside pass"  => timed(()->(@inside near[I] = active[I])),
           "one any(Bool)"          => timed(()->any(near)),
           "one fill!(Bool)"        => timed(()->fill!(near,false))]
println("N=$N on ",CUDA.name(CUDA.device()),", ",count(active)," of ",length(active)," cells active")
for (name,s) in results; @printf("%-24s %8.2f ms\n",name,1e3s); end
