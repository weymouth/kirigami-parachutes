# Free-falling parachutes, with the body motion coupled to the flow.
# Output: kirigami_N128_*_fall.jld2, used by figure_5 and figure_6 in figures/figures.jl
include(joinpath(@__DIR__,"..","src","free_fall.jl"))
using JLD2,Plots
mkpath(datadir); cd(datadir)

N = 2^7; R = 2N/3.f0; ρ = 10.f0

# single run with VTK output of the flow and body motion
let H = 0.25, θ₀ = 0.4f0
    sim = kirigami_half(N;mem=CuArray,H,fall=true,θ₀); drop!(sim,compute_parameters(N,H,θ₀,ρ;R,mem=CuArray))
    X₀(i) = sim.flow.uBC.body.X₀[i]
    motion(a) = (a.flow.f .= 0; a.flow.f[:,:,:,1] .= X₀(1); a.flow.f[:,:,:,2] .= X₀(2); a.flow.f |> Array)
    writer = vtkWriter("kirigami_N$(N)_H$(H)_fall"; attrib=Dict("u"=>vtk_u,"ω"=>vtk_ω,"λ₂"=>vtk_λ₂,"d"=>vtk_d,"motion"=>motion))
    foreach(t->(sim_step!(sim,t); save!(writer,sim)), 0.2:0.05:20.0)
    close(writer)
end

# domain sweep (the third domain, (6N,4N,3N÷2), is the H=1, θ₀=0.2 case of the sweep below)
let H = 1.f0, θ₀ = 0.2f0, times = 0.2:0.2:20.0
    params = compute_parameters(N,H,θ₀,ρ;R,mem=CuArray,T=Float32)
    for dims in ((6N,4N,3N),(8N,4N,3N))
        @show dims, θ₀, H
        sim = kirigami_half(N;mem=CuArray,H,fall=true,θ₀,dims=dims); drop!(sim,params)
        Xₘ = body_map(sim.body).x₀ # moment point in lab frame
        measure_sdf!(sim.flow.σ,sim.body,WaterLily.time(sim))
        flood(sim.flow.σ[2:end-1,2:end-1,2],clims=(-1,1))
        savefig("kirigami_N$(N)_$(dims[1])x$(dims[2])x$(dims[3])_initial.png")
        data = freefall!(sim,times)
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
        params = compute_parameters(N,H,0.2f0,ρ;R,mem=CuArray,T=Float32)
        for θ₀ in (0.f0,0.2f0,0.4f0)
            @show θ₀,H
            sim = kirigami_half(N;mem=CuArray,H,fall=true,θ₀,dims=(6N,4N,3N÷2)); drop!(sim,params)
            Xₘ = body_map(sim.body).x₀ # moment point in lab frame
            measure_sdf!(sim.flow.σ,sim.body,WaterLily.time(sim))
            flood(sim.flow.σ[2:end-1,2:end-1,2],clims=(-1,1))
            savefig("kirigami_N$(N)_H$(H)_θ$(θ₀)_initial.png")
            data = freefall!(sim,times)
            save_object("kirigami_N$(N)_H$(H)_θ$(θ₀)_fall.jld2",data)
            flood(sim.flow.u[2:end-1,2:end-1,2,1])
            scatter!([Xₘ[1]],[Xₘ[2]],markersize=5,color=:red,label=:none)
            savefig("kirigami_N$(N)_H$(H)_θ$(θ₀)_final.png")
        end
    end
end
