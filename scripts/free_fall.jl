# Free-falling parachutes, with the body motion coupled to the flow.
# Output: kirigami_N128_*_fall.jld2, used by figure_5 and figure_6 in figures/figures.jl
# NOTE: cannot run on BiotSavartBCs v1 yet; see the note on `biot_mom_step_fall!` in src/half_domain.jl
include(joinpath(@__DIR__,"..","src","half_domain.jl"))
using JLD2,Plots
mkpath(datadir); cd(datadir)

N = 2^7; R = 2N/3.f0; ρ = 10.f0

# single run with VTK output of the flow and body motion
let H = 0.25, θ₀ = 0.4f0, times = 0.2:0.05:20.0
    params = compute_paramaters!(N,H,θ₀,ρ;R,mem=CuArray,T=Float32)
    sim = kirigami_half(N;mem=CuArray,H=H,fall=true,θ₀);
    Xₘ = sim.body.a.b.map.x₀ # moment point in lab frame
    freefalling!(sim,times,params,Xₘ;save=true,name="kirigami_N$(N)_H$(H)_fall")
end

# domain sweep (the third domain, (6N,4N,3N÷2), is the H=1, θ₀=0.2 case of the sweep below)
let H = 1.f0, θ₀ = 0.2f0, times = 0.2:0.2:20.0
    params = compute_paramaters!(N,H,θ₀,ρ;R,mem=CuArray,T=Float32)
    for dims in ((6N,4N,3N),(8N,4N,3N))
        @show dims, θ₀, H
        sim = kirigami_half(N;mem=CuArray,H,fall=true,θ₀,dims=dims);
        Xₘ = sim.body.a.b.map.x₀ # moment point in lab frame
        measure_sdf!(sim.flow.σ,sim.body,WaterLily.time(sim))
        flood(sim.flow.σ[2:end-1,2:end-1,2],clims=(-1,1))
        savefig("kirigami_N$(N)_$(dims[1])x$(dims[2])x$(dims[3])_initial.png")
        data = freefalling!(sim,times,params,Xₘ)
        save_object("kirigami_N$(N)_$(dims[1])x$(dims[2])x$(dims[3])_fall.jld2",data)
        flood(sim.flow.u[2:end-1,2:end-1,2,1])
        scatter!([Xₘ[1]],[Xₘ[2]],markersize=5,color=:red,label=:none)
        savefig("kirigami_N$(N)_$(dims[1])x$(dims[2])x$(dims[3])_final.png")
    end
end

# initial angle and deployment ratio sweep
let times = 0.2:0.2:20.0
    for H in (0.25,0.5,1.0,2.0,4.0)
        # measure every time H changes (the added inertia is measured at θ₀=0.2)
        params = compute_paramaters!(N,H,0.2f0,ρ;R,mem=CuArray,T=Float32)
        for θ₀ in (0.f0,0.2f0,0.4f0)
            @show θ₀,H
            # set initial condition right
            params = Base.setindex(params, θ₀, :θ)
            sim = kirigami_half(N;mem=CuArray,H,fall=true,θ₀,dims=(6N,4N,3N÷2));
            Xₘ = sim.body.a.b.map.x₀ # moment point in lab frame
            measure_sdf!(sim.flow.σ,sim.body,WaterLily.time(sim))
            flood(sim.flow.σ[2:end-1,2:end-1,2],clims=(-1,1))
            savefig("kirigami_N$(N)_H$(H)_θ$(θ₀)_initial.png")
            data = freefalling!(sim,times,params,Xₘ)
            save_object("kirigami_N$(N)_H$(H)_θ$(θ₀)_fall.jld2",data)
            flood(sim.flow.u[2:end-1,2:end-1,2,1])
            scatter!([Xₘ[1]],[Xₘ[2]],markersize=5,color=:red,label=:none)
            savefig("kirigami_N$(N)_H$(H)_θ$(θ₀)_final.png")
        end
    end
end
