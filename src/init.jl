""" set the initial conditions of a sod shock tube """
function init_sod(g::Grid)
    mid = floor(Int64, g.xlen/2)
    rho1 = 1.0
    rho2 = 0.125
    p1 = 1.0
    p2 = 0.1

    g.rho[1:mid] .= rho1
    g.pressure[1:mid] .= p1
    g.rho[(mid+1):end] .= rho2
    g.pressure[(mid+1):end] .= p2
    g.vel .= 0.0
    return
end


""" Set the 2d initial conditions: a 2D version of the 1D sod shocktube """
function init_sod(g::Grid2d)
    mid = floor(Int64, g.xlen/2)
    rho1 = 1.0
    rho2 = 0.125
    p1 = 1.0
    p2 = 0.1

    g.rho[1:mid, :] .= rho1
    g.pressure[1:mid, :] .= p1
    g.rho[(mid+1):end, :] .= rho2
    g.pressure[(mid+1):end, :] .= p2
    g.vx .= 0.0
    g.vy .= 0.0
    return
end


""" Set the 2d sod shock tube: a ball of radius 0.2 at center with
density 1.0, pressure 1.0. Outside has density 0.125 and pressure 0.1 """
function init_ball(g::Grid2d)
    radius = 0.2
    midx = floor(Int64, g.xlen/2)
    midy = floor(Int64, g.ylen/2)
    rho1 = 1.0
    rho2 = 0.125
    p1 = 1.0
    p2 = 0.1

    cx = g.x[midx]
    cy = g.y[midy]
    for i = 1:g.xlen, j = 1:g.ylen
        if (g.x[i] - cx)^2 + (g.y[j] - cy)^2 <= radius^2
            g.rho[i, j] = rho1
            g.pressure[i, j] = p1
        else
            g.rho[i, j] = rho2
            g.pressure[i, j] = p2
        end
    end
    g.vx .= 0.
    g.vy .= 0.
    return
end


""" Kelvin-Helmholtz instability """
function init_KH(g::Grid2d)
    rho1 = 1.0
    rho2 = 2.0
    v1 = -0.5
    v2 = 0.5
    p = 2.5
    Amp = 0.01

    # g.vy .= 0.0
    g.pressure .= p
    for j = 1:g.ylen
        if abs(g.y[j] - 0.5) > 0.25
            g.rho[:, j] .= rho1
            g.vx[:, j] .= v1
            # g.u[:, j, 1] .= rho1
            # g.u[:, j, 2] .= g.u[:, j, 1] .* -0.5
        else
            g.rho[:, j] .= rho2
            g.vx[:, j] .= v2
            # g.u[:, j, 1] .= rho2
            # g.u[:, j, 2] .= g.u[:, j, 1] .* 0.5
        end
    end
    # g.u[:, :, 3] .= 0.          # vy = 0
    # pressure = ones(Float64, g.xlen, g.ylen) .* 2.5
    # g.u[:, :, 4] .= pressure ./ (g.gamma - 1)
    # # initial perturbation
    # g.u[:, :, 2] .+= 0.1 * (1 .- 2 * rand(g.xlen, g.ylen)) .* g.u[:, :, 1]
    # g.u[:, :, 3] .+= 0.1 * (1 .- 2 * rand(g.xlen, g.ylen)) .* g.u[:, :, 1]
    # initial perturbation: sin wave, Springel 2009
    σ = 0.05
    for j = 1:g.ylen
        # y = j / g.ylen
        y = g.y[j]
        # vy = 0.01 * sin.(4 * pi * collect(1:g.xlen) / g.xlen)
        # expo = - (y - 0.25)^2 / (2 * σ^2) + (y - 0.75)^2 / (2 * σ^2)
        vy = 0.1 * sin.(4 * pi * g.x) .* (exp(- (y - 0.25)^2 / (2 * σ^2)) .+ exp.(- (y - 0.75)^2 / (2 * σ^2)))
        # println("y = $y, vy max = $(maximum(vy)), expo max = $(maximum(expo))")
        # g.u[:, j, 3] .= vy .* g.u[:, j, 1]
        g.vy[:, j] .= vy
    end
    # cons2prim(g)
    return
end


""" Kelvin-Helmholtz instability """
function init_KH2(g::Grid2d)
    rho1 = 1.0
    rho2 = 2.0
    v1 = 0.5
    v2 = -0.5
    p = 2.5
    Amp = 0.01

    for j = 1:g.ylen
        if abs(g.y[j] - 0.5) > 0.25
            g.rho[:, j] .= rho1
            g.vx[:, j] .= v1
        else
            g.rho[:, j] .= rho2
            g.vx[:, j] .= v2
        end
    end
    Lx = g.x[g.xjhi] - g.x[g.xjlo]
    Ly = g.y[g.yjhi] - g.y[g.yjlo]
    for j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.vx[i, j] += Amp * (1.0 + cos(2pi * 2 * g.x[i] / Lx) * (1.0 + cos(2pi * 2 * g.y[j] / Ly)))
    end
    g.vy .= 0.0
    g.pressure .= p

    return
end


# random velocities with amplitude 0.01
function init_KH_rand(g::Grid2d)
    rho1 = 1.0
    rho2 = 2.0
    v1 = -0.5
    v2 = 0.5
    p = 2.5
    Amp = 0.01

    for j = 1:g.ylen
        if abs(g.y[j] - 0.5) > 0.25
            g.rho[:, j] .= rho1
            g.vx[:, j] .= 0.1
        else
            g.rho[:, j] .= rho2
            g.vx[:, j] .= 0.1
        end
    end
    g.vx .+= Amp .* rand(size(g.vx))
    g.vy .= Amp .* rand(size(g.vx))
    g.pressure .= p
end
