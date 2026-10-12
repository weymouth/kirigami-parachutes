include("kirigami.jl")
using CUDA

function kirigami_half(N;H=0,rings=16,U=1,a=1,Re=1e4,mem=Array,T=Float32,Ux=linear,R=T(2N/3),θ₀=zero(T),
                          dims=(3N,3N,3N÷2),ϵ=T(1/2),half_thk=ϵ+1/T(√2),fall=false,dir=1)
    H = T(H)
    x₀ = apex_offset(R,H,rings,half_thk)+R/4
    H==1/4 && (x₀ += R)
    body = kirigami_body(RigidMap(SA{T}[x₀,dims[2]/2,0],SA{T}[0,0,θ₀]),R,H,rings,half_thk)
    fall && (body = NarrowBand(body,dims;interior=false,mem))
    Ut = fall ? (0,0,0) : (i,x,t)->(i==dir ? U*Ux(a*U*t/2R) : zero(t))
    BiotSimulation(dims,Ut,R;U,ν=U*2R/Re,body,mem,T,ϵ,symmetry=(-3,))
end

drag!(sim,times,R=sim.L,x₀=SA[R,0,0];remeasure=false) = history(sim,times,4,R,x₀;remeasure)
