using Plots

function main()
    # Scalars
    rho = 1.0
    gamma = 0.0001
    u = 5.0
    v = 1.0
    N = 25
    L = 1.0
    h = L / N
    phi_L = 1.0
    phi_R = 0.0
    phi_T = 1.0
    phi_B = 0.0
    D = gamma / h           # diffusion strength
    F = rho * u             # advection strength
    Fy = rho * v            # vertical advection strength (south -> north)

    # Coordinates
    x = range(h/2, L - h/2, length=N)
    y = range(h/2, L - h/2, length=N)

    # Fields
    phi = zeros(N, N)
    phi_new = zeros(N, N)

    epsilon = 1e-6
    iterations = 0
    num_error = 1.0

    while num_error > epsilon
        num_error = 0.0

        for i in 1:N, j in 1:N
            if i >= 2 && i <= N-1 && j >= 2 && j <= N-1          # interior
                aW = D + F; aS = D + Fy; aE = D; aN = D
                aP = 4*D + F + Fy
                Su = 0.0
                phi_new[i,j] = (aW*phi[i,j-1] + aE*phi[i,j+1] + aN*phi[i-1,j] + aS*phi[i+1,j] + Su) / aP

            elseif i == 1 && j >= 2 && j <= N-1                  # top edge
                aW = D + F; aS = D + Fy; aE = D; aN = 0.0
                aP = 5*D + F + Fy
                Su = 2*D*phi_T
                phi_new[i,j] = (aW*phi[i,j-1] + aE*phi[i,j+1] + aS*phi[i+1,j] + Su) / aP

            elseif j == 1 && i >= 2 && i <= N-1                  # left edge
                aW = 0.0; aS = D + Fy; aE = D; aN = D
                aP = 5*D + F + Fy
                Su = (2*D + F)*phi_L
                phi_new[i,j] = (aE*phi[i,j+1] + aN*phi[i-1,j] + aS*phi[i+1,j] + Su) / aP

            elseif i == N && j >= 2 && j <= N-1                  # bottom edge
                aW = D + F; aS = 0.0; aE = D; aN = D
                aP = 5*D + F + Fy
                Su = (2*D + Fy)*phi_B
                phi_new[i,j] = (aW*phi[i,j-1] + aE*phi[i,j+1] + aN*phi[i-1,j] + Su) / aP

            elseif j == N && i >= 2 && i <= N-1                  # right edge
                aW = D + F; aS = D + Fy; aE = 0.0; aN = D
                aP = 5*D + F + Fy
                Su = 2*D*phi_R
                phi_new[i,j] = (aW*phi[i,j-1] + aN*phi[i-1,j] + aS*phi[i+1,j] + Su) / aP

            elseif i == 1 && j == 1                              # top left
                aW = 0.0; aS = D + Fy; aE = D; aN = 0.0
                aP = 6*D + F + Fy
                Su = 2*D*phi_T + (2*D + F)*phi_L
                phi_new[i,j] = (aE*phi[i,j+1] + aS*phi[i+1,j] + Su) / aP

            elseif i == 1 && j == N                              # top right
                aW = D + F; aS = D + Fy; aE = 0.0; aN = 0.0
                aP = 6*D + F + Fy
                Su = 2*D*phi_T + 2*D*phi_R
                phi_new[i,j] = (aW*phi[i,j-1] + aS*phi[i+1,j] + Su) / aP

            elseif i == N && j == 1                              # bottom left
                aW = 0.0; aS = 0.0; aE = D; aN = D
                aP = 6*D + F + Fy
                Su = (2*D + Fy)*phi_B + (2*D + F)*phi_L
                phi_new[i,j] = (aE*phi[i,j+1] + aN*phi[i-1,j] + Su) / aP

            else                                                 # bottom right
                aW = D + F; aS = 0.0; aE = 0.0; aN = D
                aP = 6*D + F + Fy
                Su = (2*D + Fy)*phi_B + 2*D*phi_R
                phi_new[i,j] = (aW*phi[i,j-1] + aN*phi[i-1,j] + Su) / aP
            end

            num_error += abs(phi_new[i,j] - phi[i,j])
        end

        iterations += 1
        phi .= phi_new
    end

    println("Iterations: ", iterations)

    # Plot (flip y like MATLAB's L-y)
    p = contourf(collect(x), L .- collect(y), phi, levels=128, linewidth=0,
                 colorbar=true, xlabel="x", ylabel="y")
    display(p)

    # Centerline data (original was commented out)
    # mid = (N + 1) ÷ 2
    # plot(phi[:, mid], 1 .- y)
    # plot(x, phi[mid, :])
    return phi
end

phi = main()
