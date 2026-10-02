include("half_domain.jl")

function compute_parameters(N,H,θ₀,ρ;rings=16,R=2N/3,U=1,mem=CuArray,T=Float32)
    R,H,θ₀,ρ,U = T(R),T(H),T(θ₀),T(ρ),T(U)
    sim(;kw...) = kirigami_half(N;R,T,mem,H,dims=(3N,3N,3N÷2),kw...)
    force(sim) = T.(WaterLily.pressure_force(sim))
    added_mass(dir) = (s = sim(;dir); sim_step!(s;remeasure=false); -2force(s)[dir]/outer_radius(R)^2)
    m₁₁,m₂₂ = added_mass(1),added_mass(2)

    s = sim(fall=true)
    α=one(T); ω=α*s.flow.Δt[end]; θ=θ₀+ω*s.flow.Δt[end]
    s.body = setmap(s.body;θ=SA{T}[0,0,θ],ω=SA{T}[0,0,ω])
    sim_step!(s;remeasure=true)
    Iₐ = 2T.(WaterLily.pressure_moment(body_map(s.body).x₀,s))[3]
    apply!(x->-x[1],s.flow.p)
    m = -2force(s)[1]

    (m=ρ*m, g=SA{T}[-U^2/R,0,0], mₐ=SA{T}[m₁₁*R^3,m₂₂*R^3,0], Iₘ=ρ*body_inertia(R,H,rings,ρ)/2, Iₐ)
end

struct Falling{V,A,T,P} <: Function
    U::V; a::A; t₀::T; body::P
end
(f::Falling)(i,x,t) = f.U[i] + f.a[i]*(t - f.t₀)

drop!(sim,body) = (z = zero(body.g); T = eltype(z);
    sim.flow = WaterLily.setproperties(sim.flow;uBC=Falling(z,z,zero(T),(;body...,α=zero(T),X₀=z))))

@inline @fastmath rotate(v,θ::T) where T = SA{T}[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]*v

function WaterLily.measure!(sim::Simulation,t=sum(sim.flow.Δt))
    if sim.flow.uBC isa Falling
        fall = sim.flow.uBC
        T, Δt, map = eltype(sim.flow.p), sim.flow.Δt[end], body_map(sim.body)
        (;m,mₐ,Iₘ,Iₐ,g,α,X₀) = fall.body
        vel,acc,θ,ω = -fall.U,-fall.a,map.θ[3],map.ω[3]
        force = rotate(-T.(WaterLily.total_force(sim))+m.*g,-θ)
        moment = -T.(WaterLily.pressure_moment(map.x₀,sim))[3]
        acc = rotate((force-mₐ.*rotate(acc,-θ))./(m.+mₐ),θ).*SA{T}[1,1,0]
        vel += Δt*acc; X₀ += Δt*vel
        α = (moment-α*Iₐ)/(Iₘ+Iₐ); ω += Δt*α; θ += Δt*ω
        sim.body = setmap(sim.body;θ=SA{T}[0,0,θ],ω=SA{T}[0,0,ω])
        sim.flow = WaterLily.setproperties(sim.flow;uBC=Falling(-vel,-acc,t,(;fall.body...,α,X₀)))
    end
    WaterLily.measure!(sim.flow,sim.body;t,ϵ=sim.ϵ)
    WaterLily.update!(sim.pois)
end

function freefall!(sim,times;R=sim.L)
    data = NamedTuple[]
    for t in times
        @show t; flush(stdout)
        sim_step!(sim,t)
        fall, map = sim.flow.uBC, body_map(sim.body)
        maximum(abs,fall.U) > 10sim.U && break
        push!(data,(;t,coefficients(sim,4,R,map.x₀)...,u₁=-fall.U[1],u₂=-fall.U[2],a₁=-fall.a[1],a₂=-fall.a[2],
                     θ=map.θ[3],ω=map.ω[3],fall.body.α))
    end
    Table(data)
end
