using Plots

function main()
    # Problem definition
    N = 5               # grid points
    u = 2.5
    L = 1.0
    rho = 1.0
    gamma = 0.1

    h = L / N
    x = range(h/2, L - h/2, length=N)

    F = rho * u         # advection strength
    D = gamma / h       # diffusion strength
    Pe = F / D          # Peclet number
    println("Pe = ", Pe)

    phi = zeros(N)
    phi0 = 1.0          # left boundary value
    phiL = 0.0          # right boundary value
    phi_new = copy(phi)

    epsilon = 1e-7
    it_count = 0
    residual = 1.0

    its = Int[]
    res = Float64[]

    # Main loop (upwind)
    while residual > epsilon
        for index in 2:(N-1)
            aE = D
            aW = D + F
            aP = aE + aW
            phi_new[index] = (aE*phi[index+1] + aW*phi[index-1]) / aP
        end

        # Left boundary node
        index = 1
        aE = D
        aW = 0.0
        a0 = 2*D + F
        aP = aE + aW + a0
        phi_new[index] = (aE*phi[index+1] + a0*phi0) / aP

        # Right boundary node
        index = N
        aE = 0.0
        aW = D + F
        aL = 2*D
        aP = aE + aW + aL
        phi_new[index] = (aW*phi[index-1] + aL*phiL) / aP

        it_count += 1
        residual = sum(abs.(phi_new .- phi))
        push!(its, it_count)
        push!(res, residual)

        phi .= phi_new
    end

    println("Iterations: ", it_count)

    p1 = scatter(its, res, yscale=:log10, xlabel="iteration", ylabel="residual",
                 legend=false, title="Convergence")
    p2 = plot(x, phi, linestyle=:dash, marker=:square, markersize=5, linewidth=1,
              xlabel="x", ylabel="phi", legend=false)
    display(p1)
    display(p2)
    return x, phi
end

x, phi = main()
