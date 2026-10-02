include("common.jl")

linear(t) = min(t,one(t))
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=1)

ring_spacing(R,H,rings) = (R/rings, R*H/rings^2)
apex_offset(R,H,rings,half_thk) = max(R*(1-H)/2,R/rings+half_thk-min(0,R*H))

function ring_sdf((x,y,z),R₀,R₁,x₀,x₁,ϕ,δR,half_thk)
    r,θ = hypot(y,z),atan(z,y)
    δx = x₀+tanh(π*r/δR)*(x₁-x₀)*(1+cos(4θ+ϕ))/2
    hypot(x-δx,r-clamp(r,R₀+half_thk,R₁-half_thk))-half_thk
end

function kirigami_body(map,R,H,rings,half_thk)
    δR,δH = ring_spacing(R,H,rings)
    ring(R₀,R₁,x₀,x₁,ϕ) = AutoBody((x,t)->ring_sdf(x,R₀,R₁,x₀,x₁,ϕ,δR,half_thk),map)
    H == 0 && return ring(0,R,0,0,0)
    sum(i -> ring(δR*(i-1), δR*i, δH*(i-1)^2, δH*i^2, typeof(R)(π)*(i%2)), 1:rings)
end

body_map(body::WaterLily.SetBody) = body_map(body.a)
body_map(body::AutoBody) = body.map

outer_radius(R) = R+typeof(R)(1/2)+1/typeof(R)(√2)
thickness(T) = 1+2/T(√2)

ring_mass(R₀,R₁,t,ρ) = 2ρ*π*t*(R₁-R₀)
ring_inertia(R₀,R₁,t,ρ,x₀) = (m = ring_mass(R₀,R₁,t,ρ); m*(R₁^2+R₀^2)/4 + m*x₀^2)
function body_inertia(R,H,rings,ρ)
    δR,δH = ring_spacing(R,H,rings)
    sum(i -> ring_inertia(δR*(i-1), δR*i, thickness(typeof(R)), ρ, δH*i^2-δH*(i-1)^2), 1:rings)
end

function coefficients(sim,scale,R,x₀)
    T = eltype(sim.flow.p)
    Cd,Cl = -scale*T.(WaterLily.total_force(sim))[1:2]/R^2
    Cm = scale*T.(WaterLily.pressure_moment(x₀,sim))[3]/R^3
    (;Cd,Cl,Cm)
end

history(sim,times,scale,R,x₀;remeasure=false) = map(times) do t
    @show t; flush(stdout)
    sim_step!(sim,t;remeasure)
    (;t,coefficients(sim,scale,R,x₀)...)
end |> Table
