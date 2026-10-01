# Sphere validation case, shared by scripts/sphere.jl and figures/figures.jl
include("common.jl")

function make_sphere(domain; N=2^6, R=N÷3, U=1, Re=3700, T=Float32, mem = Array)
    body = AutoBody((x,t)->√sum(abs2,x .- domain[2]÷2)-R)
    BiotSimulation(domain, (U,0,0), R; ν=U*R/Re, body, T, mem)
end
