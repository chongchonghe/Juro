include("tools/SodShockTubeExact.jl")

using .SodShockTube
using Printf
using Plots


""" make the plot or update the animation, 1D case """
function plot_curve(g::Grid; fn="t.png", is_save=true)
    # calculate rho, u, p e
    x = g.x[g.jlo:g.jhi]
    data = zeros(g.nx, 4)
    data[:, 1] .= g.rho[g.jlo:g.jhi]
    data[:, 2] .= g.pressure[g.jlo:g.jhi]
    data[:, 3] .= g.vel[g.jlo:g.jhi]
    data[:, 4] .= g.epsilon[g.jlo:g.jhi]
    # data[:, 4] .= g.E[g.jlo:g.jhi]

    # exact solutions
    # Set up a shock tube problem
    problem = ShockTubeProblem(
        geometry = (0.0, 1.0, 0.5), # left edge, right edge, initial shock location
        left_state = (ρ = 1.0, u = 0.0, p = 1.0),
        right_state = (ρ = 0.125, u = 0.0, p = 0.1),
        t = g.t,
        γ = 1.4
    )
    positions, regions, values = solve(problem, x)

    p = scatter(x, data, layout=4,
                # title=["density" "pressure" "velocity" "specific energy"],
                ms=1, legend=false, xlabel="x",
                ylabel=["rho" "p" "vel" "epsilon"], xlim=[0, 1],
                ylim=[(0., 1.1) (-0., 1.2) (-.2, 1) (1.6, 3.)],
                # ylim=[(0., 1.1) (-0., 1.2) (-.2, 1) (-0.1, 2.6)],
                dpi=300, title=@sprintf("t = %.04f", g.t))
    # annotate!(1.25, 1.25, "t = $(g.t)")

    exact_e = values.p ./ (g.gamma - 1) ./ values.ρ
    plot!(x, hcat(values.ρ, values.p, values.u, exact_e), layout=4,
          color=:blue)

    # find the square difference
    relerror = sqrt(sum((data[:, 1] .- values.ρ).^2)) / sum(values.ρ)
    println("Relative error on rho is $(relerror)")

    # # new style
    # p1 = scatter(x, g.rho[g.jlo:g.jhi], ylabel="rho", ylim=(0., 1.1), legend=false, ms=1)
    # p2 = scatter(x, g.pressure[g.jlo:g.jhi], ylabel="p", ylim=(0, 1.2), legend=false, ms=1)
    # p3 = scatter(x, g.vel[g.jlo:g.jhi], ylabel="vel", ylim=(-0.2, 1), legend=false, ms=1)
    # p4 = scatter(x, g.epsilon[g.jlo:g.jhi], ylabel="epsilon", ylim=(1.6, 3), legend=false, ms=1)
    # plot(p1, p2, p3, p4, dpi=300, layout=@layout [a b; c d])

    if is_save
        savefig(fn)
    else
        return p
    end
end


""" Plot 1D curve from a 2d simulation """
function plot_curve(g::Grid2d; fn="t.png", is_save=true)
    g1 = twod2oned(g)
    plot_curve(g1; fn=fn, is_save=is_save)
end


""" Plot heatmap of the density for a 2d simulation """
function plot_heat(g::Grid2d; fn="heat.png", is_save=true, zmin=0., zmax=2.)
    # calculate rho, u, p e
    x = g.x[g.xjlo:g.xjhi]
    y = g.y[g.yjlo:g.yjhi]
    z = g.rho[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    thesize = (800, 800)
    p0 = heatmap(y, x, z, dpi=300, size=thesize, clim=(zmin, zmax),
                 xlim=(0, 1), ylim=(0, 1), showaxis=false, c=:thermal,
                 aspectratio=:equal)
    if is_save
        savefig(fn)
    else
        return p0
    end
end


""" Plot curve for a 1d simulation and heat for a 2d simulation """
function plot_curve_or_heat(g::Grid; fn="t.png", is_save=true)
    plot_curve(g; fn=fn, is_save=is_save)
end


""" Plot curve for a 1d simulation and heat for a 2d simulation """
function plot_curve_or_heat(g::Grid2d; fn="t.png", is_save=true)
    plot_heat(g; fn=fn, is_save=is_save)
end


""" Plot heatmap of the density, 2D case """
function plot_heat_four_panels(g::Grid2d; fn="heat.png", is_save=true, zmin=0., zmax=2.)
    # calculate rho, u, p e
    x = g.x[g.xjlo:g.xjhi]
    y = g.y[g.yjlo:g.yjhi]
    z1 = g.rho[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    z2 = g.vx[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    z3 = g.pressure[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    z4 = g.vy[g.xjlo:g.xjhi, g.yjlo:g.yjhi]'
    # thesize = (1600, 1600)
    # p0 = heatmap(y, x, z1, dpi=300,
    #              layout=4, title=["rho", "vx", "pressure", "vy"],
    #              size=thesize,
    #              clim=[(zmin, zmax), (-2, 2), (0, 3), (-2, 2)],
    #              xlim=(0, 1), ylim=(0, 1), showaxis=false, c=:thermal,
    #              aspectratio=:equal)
    # new style
    p1 = heatmap(y, x, z1, title="rho", xlim=(0, 1), ylim=(0, 1), showaxis=false,
                 c=:thermal, aspectratio=:equal, clim=(0, 1.1))
    p2 = heatmap(y, x, z2, title="vx", xlim=(0, 1), ylim=(0, 1), showaxis=false,
                 c=:thermal, aspectratio=:equal, clim=(-1, 1))
    p3 = heatmap(y, x, z3, title="pressure", xlim=(0, 1), ylim=(0, 1), showaxis=false,
                 c=:thermal, aspectratio=:equal, clim=(0, 2))
    p4 = heatmap(y, x, z4, title="vy", xlim=(0, 1), ylim=(0, 1), showaxis=false,
                 c=:thermal, aspectratio=:equal, clim=(-1, 1))
    p0 = plot(p1, p2, p3, p4, size=(1800, 1800), dpi=300, layout=@layout [a b; c d])
    if is_save
        savefig(fn)
    else
        return p0
    end
end
