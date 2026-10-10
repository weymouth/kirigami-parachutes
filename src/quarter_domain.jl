include("kirigami.jl")

function kirigami_quarter(N;H=0,rings=16,U=1,Re=1e4,mem=Array,T=Float32,Ux=linear,R=T(2N/3),ϵ=T(1/2),half_thk=max(ϵ+1/T(√2),R/100),fall=false)
    H = T(H)
    x₀ = apex_offset(R,H,rings,half_thk)
    body = kirigami_body((x,t)->x-SA[x₀,0,0],R,H,rings,half_thk)
    Ut = fall ? (0,0,0) : (i,x,t)->(i==1 ? U*Ux(U*t/2R) : zero(t))
    BiotSimulation((3N,N,N),Ut,R;U,ν=U*2R/Re,body,mem,T,ϵ,symmetry=(-2,-3))
end

drag!(sim,times,R=sim.L,x₀=SA[R,0,0];remeasure=false) = history(sim,times,8,R,x₀;remeasure)
