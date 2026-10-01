# Kirigami parachute with a prescribed axial acceleration, simulated on a quarter domain with
# symmetry planes in y and z. Include either this file or half_domain.jl in a Julia session,
# not both: each one overwrites `BiotSavartBCs.symmetry` and defines its own `drag!`.
include("common.jl")

# linear acceleration profile
linear(t)=min(t,one(t))
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=1) # good idea when accelerating from rest

function kirigami_quarter(N;H=0,rings=16,U=1,Re=1e4,mem=Array,T=Float32,Ux=linear,R=T(2N/3),ϵ=T(1/2),half_thk=max(ϵ+1/T(√2),0.01R),fall=false)
    δR = R/rings; δH = R*H/rings^2; x₀ = max(R*(1-H)/2,δR+half_thk-min(0,R*H))
    @inline mapped(f) = AutoBody(f,(x,t)->x-SA[x₀,0,0])
    @inline ring(R₀,R₁,x₀,x₁,ϕ) = mapped() do (x,y,z),t
        r,θ = hypot(y,z),atan(z,y)
        δx = x₀+tanh(π*r/δR)*(x₁-x₀)*(1+cos(4θ+ϕ))/2
        hypot(x-δx,r-clamp(r,R₀+half_thk,R₁-half_thk))-half_thk
    end
    body = sum(i -> ring(δR*(i-1), δR*i, δH*(i-1)^2, δH*i^2, π*(i%2)), 1:rings)
    H == 0 && (body = ring(0,R,0,0,0))
    Ut = fall ? (0,0,0) : (i,x,t)->(i==1 ? U*Ux(U*t/2R) : zero(t)) # velocity BC
    BiotSimulation((3N,N,N),Ut,R;U,ν=U*2R/Re,body,mem,T,ϵ,nonbiotfaces=(-2,-3))
end

import BiotSavartBCs: interaction,symmetry,image
@inline function symmetry(ω,T,args...) # overwrite to add image influences
    T₂,sgn₂ = image(T,size(ω),-2)  # image target and sign in y
    T₃,sgn₃ = image(T,size(ω),-3)  # image target and sign in z
    T₂₃,_   = image(T₃,size(ω),-2) # image of image!
    # Add up the four contributions
    return interaction(ω,T,args...)+sgn₃*interaction(ω,T₃,args...)+
     sgn₂*(interaction(ω,T₂,args...)+sgn₃*interaction(ω,T₂₃,args...))
end

# force and moment coefficient histories (quarter-domain normalization)
drag!(sim,times,R=sim.L,x₀=SA[R,0,0];remeasure=false) = map(times) do t
    @show t; flush(stdout)
    sim_step!(sim,t;remeasure)
    Cd,Cl = -8WaterLily.total_force(sim)[1:2]/R^2
    Cm = 8WaterLily.pressure_moment(x₀,sim)[3]/R^3
    (;t,Cd,Cl,Cm)
end |> Table
