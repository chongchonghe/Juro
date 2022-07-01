""" Solver """

# LAX scheme
function lax(g::Grid)
    fu = calc_flux(g)
    for j = 1:3, i = g.jlo:g.jhi
        # calculate L(u)
        g.lu[i, j] = -0.5 / g.dx * (fu[i + 1, j] - fu[i - 1, j])
        # replace g.u with (g.u[i-1, j] + g.u[i+1, j]) / 2
        g.u[i, k] = 0.5 * (g.u[i-1, k] + g.u[i+1, k])
    end
    # @. g.u[g.jlo:g.jhi, :] = 0.5 * (g.u[g.jlo-1:g.jhi-1, :] + g.u[g.jlo+1:g.jhi+1, :])
    return g.lu
end


function lax(g::Grid2d)
    F, G = calc_flux(g)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        # calculate L(u)
        g.lu[i, j, k] = -0.5 / g.dx * (F[i + 1, j, k] - F[i - 1, j, k]) - 0.5 / g.dy * (G[i, j + 1, k] - G[i, j - 1, k])
        # replace g.u with (g.u[i-1, j] + g.u[i+1, j]) / 2
        g.u[i, j, k] = 0.25 * (g.u[i-1, j, k] + g.u[i+1, j, k] + g.u[i, j-1, k] + g.u[i, j+1, k])
    end
    return g.lu
end


""" HLL Riemann solver, 1st order. Piecewise constant construction is applied
where the left and right states on the interface are simply
    U_{i+1/2},L = u_{i}
    U_{i+1/2},R = u_{i+1}
The wave-speeds are estimated via the simple expression:
    S_L = min{U_L - a_L, U_R - a_R}
    S_R = max{U_L + a_L, U_R + a_R}
"""

# min(uL - aL, uR - aR) and max(uL + aL, uR + aR)
function wave_speed_estimates_1(g::Grid)
    plus = g.vel .+ g.cs
    minus = g.vel .- g.cs
    lam_plus = similar(plus)
    lam_minus = similar(plus)
    for i = g.jlo-1:g.jhi
        lam_plus[i] = max(plus[i], plus[i+1])
        lam_minus[i] = min(minus[i], minus[i+1])
    end
    return lam_minus, lam_plus
end


# uL - aL and uR + aR
function wave_speed_estimates_0(g::Grid)
    lam_minus = similar(g.cs)
    lam_plus = similar(g.cs)
    for i = g.jlo-1:g.jhi
        lam_minus[i] = g.vel[i] - g.cs[i]
        lam_plus[i] = g.vel[i+1] + g.cs[i+1]
    end
    return lam_minus, lam_plus
end


# Roe 1981 JCP 43:357–372
function wave_speed_estimates_2(g::Grid)
    # TODO: finish this.
    uhat = @. (sqrt(g.rhoL) * g.velL + sqrt(g.rhoR) * g.velR) / (sqrt(g.rhoL) + sqrt(g.rhoR))
    HL = @. (g.pressureL / (g.gamma - 1) + 0.5 * g.rhoL * g.velL^2 + g.pressureL) / g.rhoL
    HR = @. (g.pressureR / (g.gamma - 1) + 0.5 * g.rhoR * g.velR^2 + g.pressureR) / g.rhoR
    Hhat = @. (sqrt(g.rhoL) * HL + sqrt(g.rhoR) * HR) / (sqrt(g.rhoL) + sqrt(g.rhoR))
    ahat = @. sqrt(g.lambda * (Hhat - 0.5 * uhat^2))
    lam_minus = uhat - ahat
    lam_plus = uhat + ahat
    return lam_minus, lam_plus
end


function hll1st(g::Grid)
    fu = calc_flux(g)
    # calculate the eigenvalues
    alpha_plus = similar(g.vel)
    alpha_minus = similar(g.vel)
    # this is equivalent to the "if 0 <= SL, if 0 > SR" version in Toro's book
    # -------- OLD --------
    # lam_plus = g.vel .+ g.cs
    # lam_minus = g.vel .- g.cs
    # for i = g.jlo-1:g.jhi
    #     alpha_plus[i] = max(0.0, lam_plus[i], lam_plus[i+1])
    #     alpha_minus[i] = max(0.0, -lam_minus[i], -lam_minus[i+1])
    # end
    # -------- END --------
    # -------- NEW --------
    lam_minus, lam_plus = wave_speed_estimates_1(g)
    for i = g.jlo-1:g.jhi
        alpha_plus[i] = max(0.0, lam_plus[i])
        alpha_minus[i] = max(0.0, -lam_minus[i])
    end
    # -------- END --------
    # calculate the difference of flux
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = (alpha_plus[i] * fu[i, k] + alpha_minus[i] * fu[i+1, k]
                        - alpha_plus[i] * alpha_minus[i] *
                            (g.u[i+1, k] - g.u[i, k])) /
                            (alpha_plus[i] + alpha_minus[i])
    end
    for k = 1:3, i = g.jlo:g.jhi
        g.lu[i, k] = - (g.fhll[i, k] - g.fhll[i-1, k]) / g.dx
    end
    i = g.xmid
    @debug "\ng.t = $(g.t)"
    @debug "all the terms in fhll:"
    @debug alpha_plus[i], fu[i, 2], alpha_minus[i], fu[i+1, 2], \
        g.u[i+1, 2], g.u[i, 2]
    @debug "fhll[$(i), $(2)]: $(g.fhll[i, 2])"
    return g.lu
end


""" HLL, 1st order. Modified based on Toro's book and Springel's notes """
function hll1st_toro(g::Grid)
    # TODO. Start from here. See Section 10.5.1 of Toro's book

    fu = calc_flux(g)
    # calculate the eigenvalues
    lam_plus = g.vel .+ g.cs
    lam_minus = g.vel .- g.cs
    alpha_plus = similar(g.vel)
    alpha_minus = similar(g.vel)
    # this is equivalent to the "if 0 <= SL, if 0 > SR" version in Toro's book
    for i = g.jlo-1:g.jhi
        alpha_plus[i] = max(0.0, lam_plus[i], lam_plus[i+1])
        alpha_minus[i] = max(0.0, -lam_minus[i], -lam_minus[i+1])
    end
    # calculate the difference of flux
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = (alpha_plus[i] * fu[i, k] + alpha_minus[i] * fu[i+1, k]
                        - alpha_plus[i] * alpha_minus[i] *
                            (g.u[i+1, k] - g.u[i, k])) /
                            (alpha_plus[i] + alpha_minus[i])
    end
    for k = 1:3, i = g.jlo:g.jhi
        g.lu[i, k] = - (g.fhll[i, k] - g.fhll[i-1, k]) / g.dx
    end
    i = g.xmid
    @debug "\ng.t = $(g.t)"
    @debug "all the terms in fhll:"
    @debug alpha_plus[i], fu[i, 2], alpha_minus[i], fu[i+1, 2], \
        g.u[i+1, 2], g.u[i, 2]
    @debug "fhll[$(i), $(2)]: $(g.fhll[i, 2])"
    return g.lu
end


function minmod(x, y, z)
    0.25 * abs(sign(x) + sign(y)) * (sign(x) + sign(z)) *
        min(abs(x), abs(y), abs(z))
end


""" Interpolate prims (rhoL, rhoR, velL, velR, pressureL, pressureR),
which then update the cons (uL, uR) """
function interpolate(g::Grid, theta::Float64=1.5)
    cdiff = similar(g.prims)
    for k = 1:3
        c = @view g.prims[:, k]     # c = rho, vel, pressure for k = 1, 2, 3
        for i = g.jlo-1:g.jhi+1
            cdiff[i, k] = 0.5 * minmod(theta * (c[i] - c[i-1]),
                                       0.5 * (c[i+1] - c[i-1]),
                                       theta * (c[i+1] - c[i]))
        end
    end
    for k = 1:3, i = g.jlo-1:g.jhi
        g.primsL[i, k] = g.prims[i, k] + cdiff[i, k]
        g.primsR[i, k] = g.prims[i+1, k] - cdiff[i+1, k]
    end
    for i = g.jlo-1:g.jhi
        g.csL[i] = prim2cs(g.rhoL[i], g.pressureL[i], g.gamma)
        g.csR[i] = prim2cs(g.rhoR[i], g.pressureR[i], g.gamma)
    end
    # update uR and uL
    # NOTE: g.rhoL, g.velL, and g.pressureL are links to the rows in
    # g.primsL, so they update along with it.
    prim2cons!(g.rhoL, g.velL, g.pressureL, g.uL, g.gamma) # no potential sqrt error
    prim2cons!(g.rhoR, g.velR, g.pressureR, g.uR, g.gamma)
end

""" Second-order HLL. This should not change g.u """
function hll2nd(g::Grid)
    interpolate(g)
    alpha_plus = max.(0.0, g.velL .+ g.csL, g.velR .+ g.csR)
    alpha_minus = max.(0.0, -g.velL .+ g.csL, -g.velR .+ g.csR)

    # save some memory by calculating fhll step by step
    calc_flux!(g.rhoL, g.velL, g.pressureL, g.gamma, g.fu)
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = alpha_plus[i] * g.fu[i, k]
    end
    calc_flux!(g.rhoR, g.velR, g.pressureR, g.gamma, g.fu)
    for k = 1:3, i = g.jlo-1:g.jhi
        g.fhll[i, k] = (g.fhll[i, k] + alpha_minus[i] * g.fu[i, k] - alpha_plus[i] * alpha_minus[i] * (g.uR[i, k] - g.uL[i, k])) / (alpha_plus[i] + alpha_minus[i])
    end
    # calculate L(u)
    for k = 1:3, i = g.jlo:g.jhi
        g.lu[i, k] = - (g.fhll[i, k] - g.fhll[i-1, k]) / g.dx
    end

    @debug "\ng.t = $(g.t)"
    @debug "all the terms in fhll:"
    i = g.xmid
    k = 2
    @debug alpha_plus[i], g.fuL[i, k], alpha_minus[i], g.fuR[i, k], g.uR[i, k], g.uL[i, k]
    @debug "fhll[$(i), $(k)]: $(g.fhll[i, k])"

    return g.lu
end


""" HLL, 1st order. Based on Chapter 7.2.1 of Springel's class notes """
function hll1st(g::Grid2d)
    F, G = calc_flux(g)
    # vx
    lam_plus = g.vx .+ g.cs
    lam_minus = g.vx .- g.cs
    alpha_plus = similar(g.vx)
    alpha_minus = similar(g.vx)
    # this is equivalent to the "if 0 <= SL, if 0 > SR" version in Toro's book
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        alpha_plus[i, j] = max(0.0, lam_plus[i, j], lam_plus[i+1, j])
        alpha_minus[i, j] = max(0.0, -lam_minus[i, j], -lam_minus[i+1, j])
    end
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.fhll[i, j, k] = (alpha_plus[i, j] * F[i, j, k] +
            alpha_minus[i, j] * F[i+1, j, k] -
            alpha_plus[i, j] * alpha_minus[i, j] *
            (g.u[i+1, j, k] - g.u[i, j, k])) /
            (alpha_plus[i, j] + alpha_minus[i, j])
    end
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] = - (g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx
    end
    # vy
    lam_plus = g.vy .+ g.cs
    lam_minus = g.vy .- g.cs
    alpha_plus = similar(g.vy)
    alpha_minus = similar(g.vy)
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        alpha_plus[i, j] = max(0.0, lam_plus[i, j], lam_plus[i, j+1])
        alpha_minus[i, j] = max(0.0, -lam_minus[i, j], -lam_minus[i, j+1])
    end
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.fhll[i, j, k] = (alpha_plus[i, j] * G[i, j, k] +
            alpha_minus[i, j] * G[i, j+1, k] -
            alpha_plus[i, j] * alpha_minus[i, j] *
            (g.u[i, j+1, k] - g.u[i, j, k])) /
            (alpha_plus[i, j] + alpha_minus[i, j])
    end
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dy
    end
    return g.lu
end


function interpolate_x(g::Grid2d)
    theta = 1.5
    cdiff = similar(g.prims)
    for k = 1:4
        c = @view g.prims[:, :, k]     # c = rho, vx, vy, pressure for k = 1, 2, 3, 4
        for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi+1
            cdiff[i, j, k] = 0.5 * minmod(theta * (c[i, j] - c[i-1, j]), 0.5 * (c[i+1, j] - c[i-1, j]), theta * (c[i+1, j] - c[i, j]))
        end
    end
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.primsL[i, j, k] = g.prims[i, j, k] + cdiff[i, j, k]
        g.primsR[i, j, k] = g.prims[i+1, j, k] - cdiff[i+1, j, k]
    end
    for j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.csL[i, j] = prim2cs(g.rhoL[i, j], g.pressureL[i, j], g.gamma)
        g.csR[i, j] = prim2cs(g.rhoR[i, j], g.pressureR[i, j], g.gamma)
    end
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma) # no potential sqrt error
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
end


function interpolate_y(g::Grid2d)
    theta = 1.5
    cdiff = similar(g.prims)
    for k = 1:4
        c = @view g.prims[:, :, k]     # c = rho, vx, vy, pressure for k = 1, 2, 3, 4
        for j = g.yjlo-1:g.yjhi+1, i = g.xjlo:g.xjhi
            cdiff[i, j, k] = 0.5 * minmod(theta * (c[i, j] - c[i, j-1]), 0.5 * (c[i, j+1] - c[i, j-1]), theta * (c[i, j+1] - c[i, j]))
        end
    end
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.primsL[i, j, k] = g.prims[i, j, k] + cdiff[i, j, k]
        g.primsR[i, j, k] = g.prims[i, j+1, k] - cdiff[i, j+1, k]
    end
    for j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.csL[i, j] = prim2cs(g.rhoL[i, j], g.pressureL[i, j], g.gamma)
        g.csR[i, j] = prim2cs(g.rhoR[i, j], g.pressureR[i, j], g.gamma)
    end
    prim2cons!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.uL, g.gamma) # no potential sqrt error
    prim2cons!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.uR, g.gamma)
end


""" Second-order HLL. This should not change g.u """
function hll2nd(g::Grid2d)

    @assert g.ng ≥ 2

    # x component
    interpolate_x(g)
    alpha_plus = max.(0.0, g.vxL .+ g.csL, g.vxR .+ g.csR)
    alpha_minus = max.(0.0, -g.vxL .+ g.csL, -g.vxR .+ g.csR)
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.fhll[i, j, k] = alpha_plus[i, j] * g.xfu[i, j, k]
    end
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo-1:g.xjhi
        g.fhll[i, j, k] = (g.fhll[i, j, k] + alpha_minus[i, j] * g.xfu[i, j, k] - alpha_plus[i, j] * alpha_minus[i, j] * (g.uR[i, j, k] - g.uL[i, j, k])) / (alpha_plus[i, j] + alpha_minus[i, j])
    end
    # calculate L(u)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] = - (g.fhll[i, j, k] - g.fhll[i-1, j, k]) / g.dx
    end

    # y component
    interpolate_y(g)
    alpha_plus = max.(0.0, g.vxL .+ g.csL, g.vxR .+ g.csR)
    alpha_minus = max.(0.0, -g.vxL .+ g.csL, -g.vxR .+ g.csR)
    calc_flux!(g.rhoL, g.vxL, g.vyL, g.pressureL, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.fhll[i, j, k] = alpha_plus[i, j] * g.yfu[i, j, k]
    end
    calc_flux!(g.rhoR, g.vxR, g.vyR, g.pressureR, g.gamma, g.xfu, g.yfu)
    for k = 1:4, j = g.yjlo-1:g.yjhi, i = g.xjlo:g.xjhi
        g.fhll[i, j, k] = (g.fhll[i, j, k] + alpha_minus[i, j] * g.yfu[i, j, k] - alpha_plus[i, j] * alpha_minus[i, j] * (g.uR[i, j, k] - g.uL[i, j, k])) / (alpha_plus[i, j] + alpha_minus[i, j])
    end
    # calculate L(u)
    for k = 1:4, j = g.yjlo:g.yjhi, i = g.xjlo:g.xjhi
        g.lu[i, j, k] += - (g.fhll[i, j, k] - g.fhll[i, j-1, k]) / g.dx
    end

    return g.lu

end
