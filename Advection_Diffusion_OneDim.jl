using Plots

function main()
    # Problem definition
    N = 20              # grid points
    u = 2.5             # velocity
    L = 1.0             # domain length
    rho = 1.0           # density
    gamma = 0.1         # diffusion coefficient

    h = L / N
    x = range(h/2, L - h/2, length=N)

    F = rho * u         # advection strength
    D = gamma / h       # diffusion strength

    phi = zeros(N)
    phi0 = 1.0          # left boundary value
    phiL = 0.0          # right boundary value
    phi_new = copy(phi)

    epsilon = 1e-7
    it_count = 0
    residual = 1.0      # must be bigger than epsilon

    its = Int[]
    res = Float64[]

    # Main loop (central differencing)
    while residual > epsilon
        for index in 2:(N-1)
            aE = D - F/2
            aW = D + F/2
            aP = aE + aW
            phi_new[index] = (aE*phi[index+1] + aW*phi[index-1]) / aP
        end

        # Left boundary node
        index = 1
        aE = D - F/2
        aW = 0.0
        a0 = 2*D + F
        aP = aE + aW + a0
        phi_new[index] = (aE*phi[index+1] + a0*phi0) / aP

        # Right boundary node
        index = N
        aE = 0.0
        aW = D + F/2
        aL = 2*D - F
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
