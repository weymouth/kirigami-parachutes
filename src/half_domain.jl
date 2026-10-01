# Kirigami parachute moving in the x-y plane (transverse acceleration, angle of attack, rotation
# and free fall), simulated on a half domain with a symmetry plane in z. Include either this file or quarter_domain.jl in a Julia session,
# not both: each one overwrites `BiotSavartBCs.symmetry` and defines its own `drag!`.
include("common.jl")
using CUDA

WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=1) # good idea when accelerating from rest
linear(t)=min(t,one(t))

# Biot-Savart momentum step with U and acceleration prescribed
# NOTE: needs `mom_project!(...; U)`, which BiotSavartBCs v1 does not provide. To be replaced
# by a reactive `uBC` (WaterLily-jl/WaterLily.jl#320) before the free-fall runs can be made.
import WaterLily: scale_u!,conv_diff!,udf!,BDIM!,CFL,mom_project!
function biot_mom_step_fall!(sim::AbstractSimulation;udf=nothing,U,kwargs...)
    a=sim.flow; b=sim.pois
    a.u⁰ .= a.u; scale_u!(a,0); t₁ = sum(a.Δt); t₀ = t₁-a.Δt[end]
    # predictor u → u'
    conv_diff!(a.f,a.u⁰,a.σ,a.λ;ν=a.ν,perdir=a.perdir)
    udf!(a,udf,a.u⁰,t₀; kwargs...)
    BDIM!(a);
    mom_project!(a,b,1,t₁;U)
    # corrector u → u¹
    conv_diff!(a.f,a.u,a.σ,a.λ;ν=a.ν,perdir=a.perdir)
    udf!(a,udf,a.u,t₁; kwargs...)
    BDIM!(a); scale_u!(a,0.5)
    mom_project!(a,b,0.5,t₁;U)
    push!(a.Δt,CFL(a))
end

# falling body acceleration term
fall!(flow,t;acceleration) = for i ∈ 1:ndims(flow.p)
    @loop flow.f[I,i] += acceleration[i] over I ∈ CartesianIndices(flow.p)
end

function kirigami_half(N;H=0,rings=16,U=1,a=1,Re=1e4,mem=Array,T=Float32,Ux=linear,R=T(2N/3),θ₀=0.f0,
                          dims=(3N,3N,3N÷2),ϵ=T(1/2),half_thk=ϵ+1/T(√2),fall=false,dir=1)
    δR = R/rings; δH = R*H/rings^2; x₀ = max(R*(1-H)/2,δR+half_thk-min(0,R*H))+0.25R
    H==0.25 && (x₀ += R) # add some extra space for the body to fall into when H is small
    @inline mapped(f) = AutoBody(f,RigidMap(SA[x₀,dims[2]/2.f0,0],SA{T}[0,0,θ₀]))
    @inline ring(R₀,R₁,x₀,x₁,ϕ) = mapped() do (x,y,z),t
        r,θ = hypot(y,z),atan(z,y)
        δx = x₀+tanh(π*r/δR)*(x₁-x₀)*(1+cos(4θ+ϕ))/2
        hypot(x-δx,r-clamp(r,R₀+half_thk,R₁-half_thk))-half_thk
    end
    body = sum(i -> ring(δR*(i-1), δR*i, δH*(i-1)^2, δH*i^2, π*(i%2)), 1:rings)
    H == 0 && (body = ring(0,R,0,0,0))
    Ut = fall ? (0,0,0) : (i,x,t)->(i==dir ? U*Ux(a*U*t/2R) : zero(t)) # velocity BC
    BiotSimulation(dims,Ut,R;U,ν=U*2R/Re,body,mem,T,ϵ,nonbiotfaces=(-3))
end

import BiotSavartBCs: interaction,symmetry,image
@inline function symmetry(ω,T,args...) # overwrite to add image influences
    T₃,sgn₃ = image(T,size(ω),-3)  # image target and sign in z
    # Add up the two contributions
    return interaction(ω,T,args...)+sgn₃*interaction(ω,T₃,args...)
end

# force and moment coefficient histories (half-domain normalization)
drag!(sim,times,R=sim.L,x₀=SA[R,0,0];remeasure=false) = map(times) do t
    @show t; flush(stdout)
    sim_step!(sim,t;remeasure)
    Cd,Cl = -4WaterLily.total_force(sim)[1:2]/R^2
    Cm = 4WaterLily.pressure_moment(x₀,sim)[3]/R^3
    (;t,Cd,Cl,Cm)
end |> Table

# inertia
@inline Izz(R₀,R₁,m) = 1/4*m*(R₁^2+R₀^2) # moment of inertia of a ring about its diameter
@inline mass(R₀,R₁,t,ρ) = 2π*ρ*t*(R₁-R₀) # mass of a ring
@inline I₆₆(R₀,R₁,t,ρ,x₀) = Izz(R₀,R₁,mass(R₀,R₁,t,ρ)) + mass(R₀,R₁,t,ρ)*x₀^2 # moment of inertia of a ring about its center

# measures the added-mass of the body, which is needed to update the acceleration in freefalling!
function compute_paramaters!(N,H,θ₀,ρ;rings=16,R=2N/3,U=1.f0,mem=CuArray,T=Float32)
    # longitudinal added mass
    sim = kirigami_half(N;R,T,mem,H=H,fall=false,θ₀=0.0,dims=(3N,3N,3N÷2),dir=1)
    sim_step!(sim;remeasure=false)
    m₁₁ = -2WaterLily.pressure_force(sim)[1]/(R+1/2+1/T(√2))^2
    # transverse added-mass
    sim = kirigami_half(N;R,T,mem,H=H,fall=false,θ₀=0.0,dims=(3N,3N,3N÷2),dir=2)
    sim_step!(sim;remeasure=false)
    m₂₂ = -2WaterLily.pressure_force(sim)[2]/(R+1/2+1/T(√2))^2
    # rotational added-mass
    sim = kirigami_half(N;R,T,mem,H=H,fall=true,θ₀=0.f0,dims=(3N,3N,3N÷2))
    # angular acceleration is not constant, so we just measure the moment at t=0 and
    # assume it is all due to added mass (not exact but should be close for small θ₀)
    α=1.0; ω=α*sim.flow.Δt[end]; θ=θ₀+ω*sim.flow.Δt[end];
    sim.body = setmap(sim.body;θ=SA{Float32}[0,0,θ],ω=SA{Float32}[0,0,ω])
    sim_step!(sim;remeasure=true)
    Xₘ = H==0 ? sim.body.map.x₀ : sim.body.a.b.map.x₀
    m₆₆ = 2WaterLily.pressure_moment(Xₘ,sim)[3]
    c₂₆ = -2WaterLily.pressure_force(sim)[2]
    # measure the volume, mass
    apply!((x)->-x[1],sim.flow.p)
    m = -2WaterLily.pressure_force(sim)[1]
    # measure mass and moment of inertia of the body itself
    δR = R/rings; δH = R*H/rings^2;
    I = sum(i -> I₆₆(δR*(i-1), δR*i, 1+2/√2, ρ, δH*i^2-δH*(i-1)^2), 1:rings)
    # update params
    return (m=ρ*m,                                # mass of body
            g=SA{Float32}[-U^2/R,0,0],            # gravity in lab frame
            mₐ=SA{Float32}[m₁₁*R^3, m₂₂*R^3, 0],  # added mass in body frame
            Iₘ=ρ*I/2.f0, Iₐ=m₆₆,                  # added moment of inertia
            θ=θ₀, ω=0.f0, α=0.f0)
end

# helper to rotate a vector
@inline @fastmath rotate(v,θ::T) where T = SA{T}[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]*v

function freefalling!(sim,times,state,Xₘ;R=sim.L,g=state.g,X₀=zero(g),vel=zero(g),acc=zero(g),
                      θ=state.θ,ω=state.ω,α=state.α,m=state.m,Iₘ=state.Iₘ,Iₐ=state.Iₐ,
                      mₐ=state.mₐ,save=false,name="kirigami_fall")
    vtk_motion(a::AbstractSimulation) = (a.flow.f .= 0; a.flow.f[:,:,:,1] .= X₀[1];
                                         a.flow.f[:,:,:,2] .= X₀[2]; a.flow.f |> Array)
    save && (writer = vtkWriter(name; attrib=Dict("u"=>vtk_u,"ω"=>vtk_ω,"λ₂"=>vtk_λ₂,"d"=>vtk_d,"motion"=>vtk_motion)))
    data = NamedTuple[] # store data
    for t in times
        while sim_time(sim) < t
            # the step we are doing and the initial angle
            Δt = sim.flow.Δt[end]
            # compute pressure force and moment in lab frame
            force = -WaterLily.total_force(sim)
            moment = -WaterLily.pressure_moment(Xₘ,sim)[3]
            # transform to body frame
            force,acc = rotate(force+m.*g, -θ),rotate(acc, -θ)
            # update linear motion in body frame, and then back to lab frame
            acc = rotate((force - mₐ.*acc)./(m .+ mₐ), θ).*SA{Float32}[1,1,0]
            vel += Δt*acc; X₀ += Δt*vel
            # update rotation ODE
            α = (moment - α*Iₐ)/(Iₘ + Iₐ)
            ω += Δt*α; θ += Δt*ω # Verlet
            # remeasure the sim
            sim.body = setmap(sim.body;θ=SA{Float32}[0,0,θ],ω=SA{Float32}[0,0,ω]) # update rotational variables
            measure!(sim)
            biot_mom_step_fall!(sim;udf=fall!,acceleration=-acc,U=-vel)
        end
        maximum(abs, vel) > 10sim.U && break # stop if it goes out of control, probably numerical instability at that point
        save && save!(writer,sim)
        println("tU/L=",round(t,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3),
                " X₁=", round(X₀[1]/sim.L,digits=3), " θ=", round(rad2deg(θ),digits=3),
                "° u₁=", round(vel[1]/sim.U,digits=3), " a₁=", round(acc[1]/(sim.U^2/sim.L),digits=3))
        flush(stdout)
        Cd,Cl = -4WaterLily.total_force(sim)[1:2]/R^2
        Cm = 4WaterLily.pressure_moment(Xₘ,sim)[3]/R^3
        push!(data, (;t,Cd,Cl,Cm,u₁=vel[1],u₂=vel[2],a₁=acc[1],a₂=acc[2],θ,ω,α))
    end
    save && close(writer)
    return Table(data)
end
