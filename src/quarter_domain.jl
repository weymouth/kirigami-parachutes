include("kirigami.jl")

function kirigami_quarter(N;H=0,rings=16,U=1,Re=1e4,mem=Array,T=Float32,Ux=linear,R=T(2N/3),ϵ=T(1/2),half_thk=max(ϵ+1/T(√2),R/100),fall=false)
    H = T(H)
    x₀ = apex_offset(R,H,rings,half_thk)
    body = kirigami_body((x,t)->x-SA[x₀,0,0],R,H,rings,half_thk)
    Ut = fall ? (0,0,0) : (i,x,t)->(i==1 ? U*Ux(U*t/2R) : zero(t))
    BiotSimulation((3N,N,N),Ut,R;U,ν=U*2R/Re,body,mem,T,ϵ,nonbiotfaces=(-2,-3))
end

# overwrites BiotSavartBCs.symmetry: include only one of quarter_domain.jl and half_domain.jl
import BiotSavartBCs: interaction,symmetry,image
@inline function symmetry(ω,T,args...)
    T₂,sgn₂ = image(T,size(ω),-2)
    T₃,sgn₃ = image(T,size(ω),-3)
    T₂₃,_   = image(T₃,size(ω),-2)
    return interaction(ω,T,args...)+sgn₃*interaction(ω,T₃,args...)+
     sgn₂*(interaction(ω,T₂,args...)+sgn₃*interaction(ω,T₂₃,args...))
end

drag!(sim,times,R=sim.L,x₀=SA[R,0,0];remeasure=false) = history(sim,times,8,R,x₀;remeasure)
