# Flow around a sphere at Re=UR/ν=3700 on three domain sizes (validation case in the appendix).
# Output: sphere_*.jld2 (drag histories) and the flow state of the smallest domain at tU/R=53
# (sphere_320x128x128_t53.jld2), used by figure_A2 in figures/figures.jl
include(joinpath(@__DIR__,"..","src","sphere.jl"))
using CUDA,JLD2
mkpath(datadir); cd(datadir)

# size
N=2^7
params = [(5N÷2,N,N) (15N÷4,3N÷2,3N÷2) (5N,2N,2N)]
for domain in params
    # make the sim
    sim = make_sphere(domain;N=N,R=44,mem=CUDA.CuArray)
    sim_step!(sim;remeasure=false)
    time = 0.1:0.1:200
    t = map(i->string(i),domain)
    # run
    drag = map(time) do tᵢ
        sim_step!(sim,tᵢ;remeasure=false)
        @show tᵢ
        # flow state of the smallest domain at tU/R=53 for the wake rendering
        domain==params[1] && tᵢ≈53 && save!("sphere_$(t[1])x$(t[2])x$(t[3])_t53.jld2",sim)
        -WaterLily.pressure_force(sim)[1]/(0.5π*sim.L^2)
    end
    jldsave("sphere_$(t[1])x$(t[2])x$(t[3]).jld2"; p=Array(sim.flow.p),
            u=Array(sim.flow.u), time=time, drag=drag)
end
