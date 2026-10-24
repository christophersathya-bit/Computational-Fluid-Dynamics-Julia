using Plots

function main()
    N = 20
    h = 1 / (N - 1)
    x = range(0, 1, length=N)

    T = zeros(N)
    T[N] = 1.0              # boundary temperature
    T_new = copy(T)

    epsilon = 1e-7
    it_count = 0
    residual = 1.0

    its = Int[]
    res = Float64[]

    while residual > epsilon
        for index in 2:(N-1)
            T_new[index] = 0.5 * (T[index+1] + T[index-1])
        end

        it_count += 1
        residual = sum(abs.(T_new .- T))
        push!(its, it_count)
        push!(res, residual)

        T .= T_new
    end

    println("Iterations: ", it_count)

    p1 = scatter(its, res, yscale=:log10, xlabel="iteration", ylabel="residual",
                 legend=false, title="Convergence")
    display(p1)

    # p2 = plot(x, T, linestyle=:dash, marker=:square, markersize=5, xlabel="x", ylabel="T")
    # display(p2)
    return x, T
end

x, T = main()
