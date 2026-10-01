# Figures of the paper made from the simulation data in data/ (see the scripts/ folder).
# Run the whole file, or include it and call the figure_n() functions individually.
include(joinpath(@__DIR__,"..","src","sphere.jl"))
using JLD2, TypedTables, GLMakie, CairoMakie, ColorSchemes, FileIO, ReadVTK, Interpolations
CairoMakie.activate!() # GLMakie is only used for the volume rendering in render_sphere_wake

function small_time(t;Re=550, k=4√(t/Re))
    t₁ = 2.257 + k - 0.141k^2 + 0.031k^3
    t₂ = (8.996 - 41k + 143.8k^2 + 45.4k^3)*t^2
    t₃ = (20.848 - 314.08k - 1851.36k^2 - 194.8k^3)*t^4
    t₄ = (28.864 + 6.272k)*t^6
    return π/√(Re*t)*(t₁ + t₂ + t₃ + t₄)
end

using Statistics
function log_decrement(x::Vector{<:Real}; window=nothing)
    n = length(x)
    # Default window: ~1/16 of signal length
    w = isnothing(window) ? max(1, round(Int, n / 16)) : window
    # Compute rolling average (centred, with edge clamping)
    rolling_mean = [mean(x[max(1, i-w):min(n, i+w)]) for i in 1:n]
    x_centered = x .- rolling_mean
    # Find local maxima in the centred signal
    peaks = findall(i -> x_centered[i] > x_centered[i-1] && x_centered[i] > x_centered[i+1], 2:n-1) .+ 1
    # Keep only positive peaks (above local mean)
    peaks = filter(i -> x_centered[i] > 0, peaks)
    length(peaks) < 2 && (println("Need at least 2 peaks above local mean to compute log decrement"); return x_centered,0)
    decrements = [log(x_centered[peaks[i]] / x_centered[peaks[i+1]]) for i in 1:length(peaks)-1]
    return x_centered, mean(decrements)
end

# these are relative to 1 CSS px
inch = 96
pt = 4/3
cm = inch / 2.54

# domain sweep
N = 2^7; times = 0.2:0.2:20.0
ρ=10.f0; R=2N/3.f0; U=1.f0 # only values H ∈ [0,1]
θ₀ = 0.2f0; H = 1.f0

function integrate(u,v,t)
    x = zeros(length(t)); y = zeros(length(t))
    for i in 2:length(t)
        dt = t[i] - t[i-1]
        x[i] = x[i-1] + u[i-1]*dt
        y[i] = y[i-1] + v[i-1]*dt
    end
    x,y
end

function figure_3()
    H, rings = 1, 8

    # --- flow slice (N=256, θ=π/8 azimuthal cut) ---
    Nflow = 256
    R = 2Nflow / 3

    vtk = VTKFile(joinpath(datadir,"vtk_data","kirigami_N$(Nflow)_H$(H)_rings$(rings)_000000.vti"))
    pd = get_point_data(vtk)
    ω = get_data_reshaped(pd["ω"])   # (3, nx, ny, nz)
    d = get_data_reshaped(pd["d"])   # (nx, ny, nz) signed distance to body

    nx, ny, nz = size(ω, 2), size(ω, 3), size(ω, 4)

    # index 1 is an unpopulated ghost layer (identically zero) along every axis - drop it
    xs, ys, zs = 2:nx, 2:ny, 2:nz
    ωy_itp = linear_interpolation((xs, ys, zs), ω[2, xs, ys, zs])
    ωz_itp = linear_interpolation((xs, ys, zs), ω[3, xs, ys, zs])
    d_itp = linear_interpolation((xs, ys, zs), d[xs, ys, zs])

    # θ=π/8: where cos(4θ+ϕ) = 0 for both ring phases (ϕ=0,π), so every ring sits at its
    # neutral radial offset instead of alternating extremes - avoids the θ=0 ring-fusing artifact
    θ = π / 8
    sθ, cθ = sincos(θ)

    rs = 2:ny   # radial samples, index units, r=0 at rs[1]
    X = [(xi - 2) / R for xi in xs]
    Rr = [(ri - 2) / R for ri in rs]

    ωθ = similar(X, length(xs), length(rs))
    body = similar(ωθ)
    for (ix, xi) in enumerate(xs), (ir, ri) in enumerate(rs)
        yc, zc = 2 + (ri - 2) * cθ, 2 + (ri - 2) * sθ
        ωθ[ix, ir] = -sθ * ωy_itp(xi, yc, zc) + cθ * ωz_itp(xi, yc, zc)
        body[ix, ir] = d_itp(xi, yc, zc)
    end

    clim = 0.2
    thresh = 0.02  # mask near-zero band: mostly far-field solver noise, not real structure
    flowdata = clamp.(ωθ, -clim, clim)
    flowdata[abs.(ωθ).<thresh] .= NaN

    # --- convergence sweep (N=64..384) ---
    Ns = (64, 96, 128, 192, 256, 384)
    data = Dict(N => load_object(joinpath(datadir,"kirigami_N$(N)_H$(H)_rings$(rings)_hist.jld2")) for N in Ns)
    Nref = last(Ns)
    Cd_ref = data[Nref].Cd
    coarse = collect(Ns[1:end-1])
    err = [sqrt(sum((data[N].Cd .- Cd_ref) .^ 2) / sum(Cd_ref .^ 2)) for N in coarse]

    # --- figure: flow on top, Cd traces + convergence below ---
    f = Figure(size=(22cm, 14cm), figure_padding=6, fontsize=9pt)

    axflow = Axis(f[2, 1:4], aspect=DataAspect(), xlabel="X/R", ylabel="r/R (θ=π/8)")
    co = contourf!(axflow, X, Rr, flowdata, levels=range(-clim, clim, 17), colormap=:RdBu,
        extendlow=:auto, extendhigh=:auto)
    contour!(axflow, X, Rr, body, levels=[0], color=:black, linewidth=1.5)
    Colorbar(f[1, 1:4], co, vertical=false, label="Azimuthal vorticity ω_θ", labelsize=12, width=Relative(0.6))
    xlims!(axflow, extrema(X)...)
    ylims!(axflow, extrema(Rr)...)

    # columns 1 and 4 are narrow spacers so the bottom pair's combined width matches axflow above
    ax1 = Axis(f[3, 2], xlabel="Convective time", ylabel="Drag coefficient")
    colors = get(ColorSchemes.Blues, range(0.35, 1.0, length=length(Ns)))
    for (c, N) in zip(colors, Ns)
        lines!(ax1, data[N].t, data[N].Cd, color=c, linewidth=2, label="R=$(round(Int, 2N/3))")
    end
    axislegend(ax1, position=:lt, rowgap=-6, patchlabelgap=2, patchsize=(16, 20),
        padding=(10, 10, -2, -2), framevisible=false)
    xlims!(ax1, 0, 3)

    Rcoarse = round.(Int, 2 .* coarse ./ 3)
    ax2 = Axis(f[3, 3], xlabel="Grid resolution R", ylabel="log₁₀(L2 error)",
        xscale=log10, xticks=(Rcoarse, string.(Rcoarse)), yticks=[-1, -1.5, -2])
    scatter!(ax2, Rcoarse, log10.(err), color=:black, markersize=10)
    lines!(ax2, Rcoarse, log10.(err), color=:black, linewidth=1)
    Nref_line = range(Rcoarse[1], Rcoarse[end], 50)
    o1 = log10.(err[1] .* (Rcoarse[1] ./ Nref_line) .^ 1)
    o2 = log10.(err[1] .* (Rcoarse[1] ./ Nref_line) .^ 2)
    lines!(ax2, Nref_line, o1, color=:gray, linestyle=:dot, linewidth=1, label="1st order")
    lines!(ax2, Nref_line, o2, color=:gray, linestyle=:dash, linewidth=1, label="2nd order")
    axislegend(ax2, position=:lb, rowgap=-6, patchlabelgap=2, patchsize=(16, 20),
        padding=(10, 10, -2, -2), framevisible=false)
    ylims!(ax2, minimum(log10.(err)) - 0.1, maximum(log10.(err)) + 0.05)

    rowgap!(f.layout, 10)
    colgap!(f.layout, 10)
    colgap!(f.layout, 1, 0)  # no gap between left spacer and ax1
    colgap!(f.layout, 3, 0)  # no gap between ax2 and right spacer
    rowsize!(f.layout, 2, Auto(1.0))  # bigger flow plot
    rowsize!(f.layout, 3, Auto(0.625))  # bottom two plots: 50% -> 62.5% (+25%) of their original size
    colsize!(f.layout, 1, Fixed(0.34cm))  # left spacer: aligns ax1 with axflow's left edge
    colsize!(f.layout, 4, Fixed(1.25cm))  # right spacer: aligns ax2 with axflow's right edge

    save(joinpath(figdir,"kirigami_convergence.png"), f, px_per_unit=300/inch)
    return f
end


function figure_4()
    inch, pt = 96, 4/3
    cm = inch / 2.54
    # added-mass and inertia data
    N = 2^8
    R=2*2^7/3.f0
    params_all = load_object(joinpath(datadir,"kirigami_parameters.jld2"))
    H = [0.0,0.25,0.5,1.0,2.0,4.0]
    mA = getindex.(params_all,:mₐ)./R^3
    m_star = getindex.(params_all,:m)./R^3
    m11 = getindex.(mA,1)
    m22 = getindex.(mA,2)
    Ia = getindex.(params_all,:Iₐ)./R^5
    Im = getindex.(params_all,:Iₘ)./R^5
    # color scheme for the rings sweep, H=1
    colors = get(ColorSchemes.amp, range(0.0, 1.0, length=6))
    # figure  and axis
    f = Figure(size=(16.4cm,5.4cm), figure_padding=6, fontsize=8pt)
    ax1 = Axis(f[1, 1], xlabel="Time (tU/R)", ylabel="Drag Coefficient")
    # ax2 = Axis(f[1, 2], xlabel="Time (tU/R)", ylabel="Drag Coefficient")
    ax2 = Axis(f[1, 3], xlabel="Deployment (H)", ylabel="Added-mass Coefficient")
    ax21 = Axis(f[1, 4], xlabel="Deployment (H)", ylabel="Lift and Moment Slope")
    ax3 = Axis(f[1, 2], xlabel="Rings (Nᵣ)", ylabel="Added-mass Coefficient")
    # analytical added-mass for N-ring disk
    Ns = 2 .^ range(0,6,200)
    lines!(ax3, Ns, π^2.0./4Ns, color=:black, linestyle=:dot, linewidth=1)
    # get results
    cs = get(ColorSchemes.hot, range(0.0, 1.0, length=9))[2:end-1]
    for (c,ms,H) in zip(cs[2:6],[:+,:circle,:rect,:diamond],(0.5,1,2,4))
        for (color,rings) = zip(colors[2:end], 4:4:20)
            data = load_object(joinpath(datadir,"kirigami_N$(N)_H$(H)_rings$(rings)_hist.jld2"))
            # H==1 && lines!(ax2,data.t,data.Cd,label="Nᵣ=$rings",linewidth=2,color=color)
            if (π^2/4rings<data.Cd[1])
                rings!=16 && scatter!(ax3, rings, data.Cd[1]; color=c, marker=ms, markersize=10)
                rings==16 && scatter!(ax3, rings, data.Cd[1]; color=c, marker=ms, markersize=10, label="H=$(H)")
            end
        end
    end
    ρ = 10
    # added mass and added inertia
    lines!(ax2,H,2m11,label="C₁₁",color=:olive)
    lines!(ax2,H,2m22,label="C₂₂",color=:maroon)
    scatter!(ax2,[0,0],[8/3,0],color=[:olive, :maroon])
    lines!(ax2,H,2Ia,label="C₆₆",color=:teal)
    scatter!(ax2,[0],[16π/45],color=:teal)
    # lines!(ax2,H,Im,label="Iₘ",color=:black, linestyle=:dash)

    lines=[]; cm_mean=[]; cl_mean=[]; cd_mean=[]
    for H in (0.0,0.25,0.5,1.0,2.0,4.0)
        data = load_object(joinpath(datadir,"kirigami_N128_H$(H)_AoA_fall.jld2"))
        # averaging window
        idx = 28 .< data.t .< 30
        # compute average
        push!(cm_mean, sum(data.Cm[idx])/length(data.Cm[idx]))
        push!(cd_mean, sum(data.Cd[idx])/length(data.Cd[idx]))
        push!(cl_mean, sum(data.Cl[idx])/length(data.Cl[idx]))
    end
    Hs = [0.0,0.25,0.5,1.0,2.0,4.0]

    println("H:    ", Hs)
    println("m*:   ", round.(2m_star.+2m22, digits=3))
    println("C_11:   ", round.(2m11, digits=3))
    println("C_26: ")
    println("C_66: ", round.(2Ia, digits=3))
    println("CL': ", round.((cl_mean./0.1), digits=3))
    println("Cm': ", round.((cm_mean./0.1), digits=3))
    println("Cd: ", round.(cd_mean, digits=3))

    scatterlines!(ax21, Hs, (cl_mean./0.1), marker=:rect, label="Cₗ'")
    scatterlines!(ax21, Hs, (cm_mean./0.1), marker=:circle, label="Cₘ'")
    # scatterlines!(ax21, Hs, cd_mean, marker=:diamond, label="Cd")
    Hi = 0:0.01:4.0; Nr=16
    lines!(ax2, Hi, π^2/16Nr .+ 0.065.*Hi.^2, color=:black, linestyle=:dot, linewidth=1)

    axislegend(ax21, position=:ct, rowgap=-6, patchlabelgap=2, patchsize=(16,20), padding=(10, 10, -2, -2))
    axislegend(ax2, position=:rt, rowgap=-6, patchlabelgap=2, patchsize=(16,20), padding=(10, 10, -2, -2))
    axislegend(ax3, position=:rt, rowgap=-6, patchlabelgap=2, patchsize=(16,20), padding=(10, 10, -2, -2))
    # xlims!(ax2,0,3); ylims!(ax2,0,8);
    xlims!(ax3,1,32); ylims!(ax3,0,2)
    # load results for H sweep
    data = load_object(joinpath(datadir,"kirigami_N$(N)_H0.0_hist.jld2"))
    lines!(ax1,data.t,data.Cd,label="H=0";color=:black)
    hlines!(ax1,[8/3],color=:black,linestyle=:dot,linewidth=1); # label="disk Ma",
    hlines!(ax1,[π^2/4/16],color=:black,linestyle=:dash,linewidth=1)# label="thin ring Ma",
    # deployment and colorscheme for that
    Hs = 2.0.^(-2:2)
    colors = get(ColorSchemes.hot, range(0.0, 1.0, length=length(Hs)+4))[2:end-1]
    for (color,H) ∈ zip(colors,Hs)
        data = load_object(joinpath(datadir,"kirigami_N$(N)_H$(H)_hist.jld2"))
        lines!(ax1,data.t,data.Cd,label="H=$(H)";linewidth=2,color=color)
        data = load_object(joinpath(datadir,"kirigami_N$(N)_H-$(H)_hist.jld2"))
        lines!(ax1,data.t,data.Cd;linewidth=2,color=color,linestyle=:dash)
    end
    xlims!(ax1,0,3.0); ylims!(ax1,-1,11)
    # colsize!(f.layout, 2, Auto(0.62))
    # colsize!(f.layout, 3, Auto(0.62))
    # colsize!(f.layout, 4, Auto(0.62))
    rowgap!(f.layout, 0)
    colgap!(f.layout, 10)
    # axislegend(ax1, position=:lt, rowgap=-6, patchlabelgap=2, patchsize=(16,20), nbanks=2, padding=(10, 10, -2, -2))
    # Label(f[1, 1, TopLeft()], "a", fontsize = 16, font = :bold, halign = :right)
    # Label(f[1, 2, TopLeft()], "b", fontsize = 16, font = :bold, halign = :right)
    # Label(f[1, 3, TopLeft()], "c", fontsize = 16, font = :bold, halign = :right)
    save(joinpath(figdir,"kirigami_Cd_time.png"), f, px_per_unit = 300/inch); f
end

function figure_5()
    f = Figure(size=(16.4cm,5.4cm), figure_padding=6, fontsize=9pt)
    ax0 = Axis(f[1:2, 1], xlabel="Y/R", ylabel="X/R", xticks=([0.0,0.3,.6],["0","0.3","0.6"]))
    # ax_inset = Axis(f[1:2, 1],width=Relative(0.35),height=Relative(0.35),halign=0.3,valign=0.2)
    ax1 = Axis(f[2, 2], xlabel="Time (tU/R)", ylabel="Coefficients")
    ax2 = Axis(f[2, 3], xlabel="Time (tU/R)", ylabel="Kinematics")
    for (dims, style) in zip(((6N,4N,3N÷2),(6N,4N,3N),(8N,4N,3N)), (:dot,:dash,:solid))
        @show dims, θ₀, H
        if dims == (6N,4N,3N÷2)
            data = load_object(joinpath(datadir,"kirigami_N$(N)_H1.0_θ0.2_fall.jld2"))
        else
            data = load_object(joinpath(datadir,"kirigami_N$(N)_$(dims[1])x$(dims[2])x$(dims[3])_fall.jld2"))
        end
        x,y = integrate(data.u₁, data.u₂, data.t)
        lines!(ax0,y,x,alpha=0.5,linewidth=2.0,linestyle=style,color=:black)
        lines!(ax1,data.t,data.Cd,alpha=0.5,linewidth=2.0,linestyle=style,color=:teal)
        lines!(ax1,data.t,data.Cm,alpha=0.5,linewidth=2.0,linestyle=style,color=:maroon)
        lines!(ax2,data.t,data.u₁,alpha=0.5,linewidth=2.0,linestyle=style,color=:firebrick)
        lines!(ax2,data.t,data.u₂,alpha=0.5,linewidth=2.0,linestyle=style,color=:olive)
        lines!(ax2,data.t,data.θ,alpha=0.5,linewidth=2.0,linestyle=style,color=:royalblue1)
    end
    xlims!(ax0,0,0.6); ylims!(ax0,-22,0)
    xlims!(ax1, 0, 20)
    xlims!(ax2,0, 20)
    rowgap!(f.layout, 0)
    colgap!(f.layout, 10)
    l1=lines!(ax1,[0],[0],linewidth=2.0,color=:black,linestyle=:solid)
    l2=lines!(ax1,[0],[0],linewidth=2.0,color=:black,linestyle=:dash)
    l3=lines!(ax1,[0],[0],linewidth=2.0,color=:black,linestyle=:dot)
    lines!(ax1,[0],[0],label="Cₓ",linewidth=2.0,color=:teal)
    lines!(ax1,[0],[0],label="Cₘ",linewidth=2.0,color=:maroon)
    lines!(ax2,[0],[0],label="u₁/U",linewidth=2.0,color=:firebrick)
    lines!(ax2,[0],[0],label="u₂/U",linewidth=2.0,color=:olive)
    lines!(ax2,[0],[0],label="θ",linewidth=2.0,color=:royalblue1)
    axislegend(ax1, position=:rc, orientation=:vertical, rowgap=-6, patchlabelgap=2, patchsize=(16,20), padding=(10, 10, 2, 2))
    axislegend(ax2, position=:rc, orientation=:vertical, rowgap=-6, patchlabelgap=2, patchsize=(16,20), padding=(10, 10, 2, 2))
    colsize!(f.layout, 1, Auto(0.5))
    Legend(f[1, 2:end], [l3, l2, l1], ["($(6N)x$(4N)x$(3N÷2))", "($(6N)x$(4N)x$(3N))", "($(8N)x$(4N)x$(3N))"],
           orientation=:horizontal, framevisible=false, patchlabelgap=2, patchsize=(26,0), padding=(10, 10, 10, 2))
    resize_to_layout!(f)
    save(joinpath(figdir,"kirigami_domain_sweep.png"), f, px_per_unit=300/inch); f
end

function figure_6()
    f = Figure(size=(16.4cm,5.4cm), figure_padding=6, fontsize=9pt)
    ax0 = Axis(f[1:2, 1], xlabel="Y/R", ylabel="X/R")
    ax1 = Axis(f[2, 2], xlabel="Time (tU/R)", ylabel="Fall Velocity (u₁/U)")
    ax2 = Axis(f[2, 3], xlabel="Deployment (H)", ylabel="ζ(u₁/U), ū₁/U",
            xticks=([0.25,1,2,4],["0.25","1","2","4"]))
    l_H=[]; l_θ=[]
    # theta and H sweep
    for θ₀ in (0.4f0,0.2f0)
        color = ifelse(θ₀ == 0.4f0, :teal, ifelse(θ₀ == 0.2f0, :maroon, :olive))
        v_final = []; δ_final = []
        for H in (0.25,0.5,1.0,2.0,4.0)
            @show θ₀,H
            data = load_object(joinpath(datadir,"kirigami_N$(N)_H$(H)_θ$(θ₀)_fall.jld2"))
            # trajectory;
            x,y = integrate(data.u₁, data.u₂, data.t)
            lines!(ax0,y,x,alpha=0.5,linewidth=min(3,H),color=color)
            H==4.0 && (l=lines!(ax1,[0],[0],linewidth=min(H,3),color=color); push!(l_θ,l))
            # final velocity extracted from the last 10% of the signal
            push!(v_final, mean(data.u₁[end-round(Int, length(data.u₁)/10):end]))
            # v_final = sum(data.u₁[end-10:end])/10
            # scatter!(ax2, [H], [v_final], marker=:rect, color=color, markersize=8+12θ₀)
            lines!(ax1,vcat(0,data.t),vcat(0,data.u₁),linewidth=min(H,3),color=color)
            H==4.0 && (l=lines!(ax1,[0],[0],linewidth=min(H,3),color=color))
            # log decreement
            xc, δ = log_decrement(vcat(0,data.u₁); window=round(Int, length(data.u₁)/32))
            push!(δ_final, δ)
            # scatter!(ax2, [H], [δ], color=color, markersize=8+12θ₀)
        end
        scatterlines!(ax2, [0.25,0.5,1.0,2.0,4.0], v_final, color=color)
        scatterlines!(ax2, [0.25,0.5,1.0,2.0,4.0], 5δ_final./sqrt.(4π^2 .+ δ_final.^2), marker=:rect, color=color)
    end
    lines!(ax1,[0],[0],label="u₁",linewidth=2.0,linestyle=:solid,color=:black)
    lines!(ax1,[0],[0],label="u₂",linewidth=2.0,linestyle=:dash,color=:black)
    xlims!(ax0,-1.2,1.8); ylims!(ax0,-25,0)
    colsize!(f.layout, 1, Auto(10/28)) # make it 1:1 aspect ratio
    rowgap!(f.layout, 0)
    colgap!(f.layout, 10)
    xlims!(ax1,0,20); ylims!(ax1,-1.5,0.0)
    # xlims!(ax2,0,20); ylims!(ax2,-1,4)
    # xlims!(ax2,1e-1,1e1); ylims!(ax2,1e-3,1e2)
    xlims!(ax2,0,4.1); ylims!(ax2,-1.5,0.6)
    for H in [0.25,0.5,1,2,4]
        l=lines!(ax1,[0],[0],linewidth=min(H,3),color=:black); push!(l_H,l)
    end
    for (m,l) in zip([:circle, :rect], ["ū₁/U", "5ζ(u₁/U)"])
        scatterlines!(ax2,[-10],[0],linewidth=2.0,color=:black,marker=m,label=l)
    end
    # @show l_H, l_θ
    axislegend(ax2, position=:rc, orientation=:horizontal, rowgap=-6, patchlabelgap=2, patchsize=(16,20), padding=(10, 10, -2, -2))
    Legend(f[1, 3], reverse(l_θ), map(θ->"θ₀=$(θ)",[0.2,0.4]), colgap=6,
            orientation=:horizontal, framevisible=false, patchlabelgap=2, patchsize=(10,20))
    Legend(f[1, 2], l_H, map(H->"H=$(H)",[0.25,0.5,1,2,4]), colgap=6,
            orientation=:horizontal, framevisible=false, patchlabelgap=2, patchsize=(10,20))
    # Legend(f[1, 3], lines[8:end], , colgap=6,
            # orientation=:horizontal, framevisible=false, patchlabelgap=2, patchsize=(10,20))
    resize_to_layout!(f)
    save(joinpath(figdir,"kirigami_results.png"), f, px_per_unit=300/inch); f
end

function figure_A1()
    data = jldopen(joinpath(datadir,"ImpCircle_results.jld2"))["single_stored_object"]
    ω = data["ω"] # vorticity on the smallest domain at tU/D=4
    biot, refl, Gillis, Koumoutsakos = data["biot"], data["refl"], data["Gillis"], data["koumoutsakos"]
    f = Figure(size=(1000,400),figure_padding=6)
    ax0 = Axis(f[2,1],aspect = DataAspect())
    ax1 = Axis(f[2,2],ylabel="Drag coefficient",xlabel="convective time")
    l1=scatter!(ax1,koumoutsakos[1,:],koumoutsakos[2,:],marker=:rect,markersize=16,
                color=:transparent,strokecolor=:black,strokewidth=1.5)
    l2=scatter!(ax1,Gillis[1,:],Gillis[2,:],marker=:circle,markersize=16,
                color=:transparent,strokecolor=:black,strokewidth=1.5)
    l3=lines!(ax1,collect(0:0.01:0.35),small_time.(0:0.01:0.35),linewidth=2,linestyle=:dash,color=:black)
    rmap = get(ColorSchemes.Reds, range(0.0, 1.0, length=6))
    D,m_2,m_1 = 128,(4,5,6),(5,10,20,40)
    for (i,m,dat) in zip(1:4,m_1,refl)
        m2 = m/2
        mod(m2,1)==0 && (m2 = Int(m2))
        lines!(ax1,dat.t,dat.Cd,color=rmap[i+1],linewidth=2,linestyle=:dash,label="Reflection, D/W=1/$m2")
    end
    bmap = get(ColorSchemes.Blues, range(0.0, 1.0, length=5))
    for (m,dat) in zip(m_2,biot)
        m2 = m/2
        mod(m2,1)==0 && (m2 = Int(m2))
        lines!(ax1,dat.t,dat.Cd,color=bmap[m-2],linewidth=2,label="Present, D/W=1/$m2")
    end
    axislegend(ax1,position=:rb, orientation=:horizontal, nbanks=4)
    ylims!(ax1,0,1.6); xlims!(ax1,0,6)
    # flow plot
    co = contourf!(ax0,clamp.(ω[inside(ω)],-5.99,5.99),levels=union(-6:-1,1:6),colormap=:RdBu)
    Colorbar(f[1,1], co, vertical=false, label="Vorticity (normalized)",labelsize=16,width=500)
    xlims!(ax0,extrema(axes(ω[inside(ω)],1))...)
    ylims!(ax0,extrema(axes(ω[inside(ω)],2))...)
    hidespines!(ax0)
    hidedecorations!(ax0)
    colsize!(f.layout, 1, Auto(1.5))
    rowsize!(f.layout, 1, Auto(0.2))
    Legend(f[1,2], [l1,l2,l3], ["Koumoutsakos et al.","Gillis et al.","Theoretical curve"], orientation=:horizontal)
    save(joinpath(figdir,"ImpCircle_results.png"),f); f
end


# volume rendering of the sphere wake, |ω|R/U on the smallest domain at tU/R=53, with WaterLily's viz!
function ω_mag!(arr,sim)
    a = sim.flow.σ
    @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*sim.L/sim.U
    copyto!(arr, a[inside(a)])
end
function render_sphere_wake(file=joinpath(figdir,"sphere3_zoom.png"))
    N = 2^7; sim = make_sphere((5N÷2,N,N);N,R=44)
    load!(sim; fname="sphere_320x128x128_t53.jld2", dir=datadir)
    GLMakie.activate!() # volume rendering needs GLMakie
    try
        viz!(sim; f=ω_mag!, img=file, fig=Figure(size=(2926,1154)), colormap=:Blues,
             colorrange=(0,20), body_color=:white, hidedecorations=true, azimuth=1.5π, elevation=0)
    finally
        CairoMakie.activate!()
    end
    return file
end

function figure_A2()
    # load the files
    small = jldopen(joinpath(datadir,"sphere_320x128x128.jld2"))
    medium = jldopen(joinpath(datadir,"sphere_480x192x192.jld2"))
    large = jldopen(joinpath(datadir,"sphere_640x256x256.jld2"))
    # get a colorscheme
    blues = get(ColorSchemes.Blues, range(0.0, 1.0, length=8))[3:end]
    img = load(render_sphere_wake())
    f = Figure(size=(1000,300), figure_padding=5)
    drag = Axis(f[2, 1], xlabel="Convective time", ylabel="Drag coefficient")
    labels = ["3.6Dx1.5Dx1.5D","5.5Dx2.2Dx2.2D","7.2Dx2.9Dx2.9D"]
    lines = []; labs = []
    for (i,case,D) in zip([2,4,6],[small,medium,large],[128,192,256])
        t = case["time"]; idx = t .> 100
        lines!(drag, t, case["drag"], color=blues[i],alpha=0.6)
        fx, t = case["drag"][idx], t[idx]
        CD_mean = sum(fx[2:end].*diff(t))/sum(diff(t))
        println("▷ ΔT [CTU] = $(t[end]-t[1])")
        println("▷ CD mean = $CD_mean")
        l=hlines!(drag, [CD_mean], color=blues[i], linewidth=2)
        push!(lines, l)
        push!(labs, (π*44^2)/D^2)
    end
    # image is of size
    nx,ny=size(img')
    # domain is 3.6Dx1.5D, this means
    Flow = Axis(f[1:2, 2], aspect = DataAspect(), xlabel="X/R", ylabel="Y/R", xticks=(0:Int(nx÷3.6):nx,["0","2","4","6"]),
                yticks=(0:ny÷2:ny,["-1.5","0","1.5"]))
    image!(Flow, img')
    hlines!(drag, [0.394], linestyle=:dash, color=:black, label="Rodriguez et al. (DNS)")
    hlines!(drag, [0.355], linestyle=:dot, color=:black, label="Yun et al. (LES)")
    xlims!(drag,0,200); ylims!(drag,0.3,0.5)
    axislegend(drag, position=:rt, labelsize=12, rowgap=0)
    Legend(f[1,1], lines, map(i->"=$(round(labs[i];digits=3))",1:3), "Blockage ratio πR²/A", labelsize=12, nbanks=3,
            framevisible=false, position=:rt, colgap=5, titlegap=0)
    # hidespines!(ax.axis); hidedecorations!(ax.axis)
    colsize!(f.layout, 1, Relative(1/4))
    rowsize!(f.layout, 1, Auto(0.15))
    rowgap!(f.layout, 1, Relative(0.02))
    save(joinpath(figdir,"validation_sphere.png"), f); f
end


# main: make all figures when run as a script
if abspath(PROGRAM_FILE) == @__FILE__
    figure_3()  # kirigami_convergence.png
    figure_4()  # kirigami_Cd_time.png
    figure_5()  # kirigami_domain_sweep.png
    figure_6()  # kirigami_results.png
    figure_A1() # ImpCircle_results.png
    figure_A2() # validation_sphere.png
end
