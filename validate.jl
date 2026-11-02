# validate.jl — checks every converted script against independent references.
# Run:  julia --project=. validate.jl     (exit code 1 if any check fails)

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")   # headless plotting
using LinearAlgebra, Printf

const DIR = @__DIR__
const results = Tuple{String,Bool}[]

function check(name, ok, detail="")
    push!(results, (name, ok))
    @printf("[%s] %s  %s\n", ok ? "PASS" : "FAIL", name, detail)
end

# ---------- helpers ----------

# Run a solver script in a throwaway module, stdout suppressed; returns the module
function load(file)
    m = Module(:Sandbox)
    path = tempname()
    open(path, "w") do io
        redirect_stdout(io) do
            Base.include(m, joinpath(DIR, file))
        end
    end
    return m
end

# Run a script as a subprocess with canned stdin; returns stdout lines
function run_script(file; input="")
    cmd = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) $(joinpath(DIR, file))`
    out = read(pipeline(cmd; stdin=IOBuffer(input)), String)
    return split(chomp(out), '\n')
end

safely(f, name) = try
    f()
catch e
    check(name, false, "threw: " * sprint(showerror, e))
end

f(x)  = -4x^3 + 7x^2 - 3x + 9
parse_show(lines, key) = parse(Float64, split(last(filter(l -> startswith(l, key * " = "), lines)), " = ")[2])

# ---------- 1D advection-diffusion: direct solve of the same discrete equations ----------

function direct_1d(N, u, L, rho, gamma, scheme)
    h = L / N; F = rho * u; D = gamma / h
    aE, aW, a0, aL = scheme == :central ?
        (D - F/2, D + F/2, 2D + F, 2D - F) :
        (D,       D + F,   2D + F, 2D)
    A = zeros(N, N); b = zeros(N)
    for i in 1:N
        if i == 1
            A[i,i] = aE + a0;  A[i,i+1] = -aE;  b[i] = a0 * 1.0   # phi0 = 1
        elseif i == N
            A[i,i] = aW + aL;  A[i,i-1] = -aW;  b[i] = 0.0        # phiL = 0
        else
            A[i,i] = aE + aW;  A[i,i+1] = -aE;  A[i,i-1] = -aW
        end
    end
    return A \ b
end

safely("Advection_Diffusion_OneDim") do
    m = load("Advection_Diffusion_OneDim.jl")
    ref = direct_1d(20, 2.5, 1.0, 1.0, 0.1, :central)
    err = maximum(abs.(m.phi .- ref))
    check("1D central vs direct solve", err < 1e-4, @sprintf("max err = %.2e", err))

    x = collect(m.x); Pe = 2.5 * 1.0 / 0.1
    exact = 1 .- (exp.(Pe .* x) .- 1) ./ (exp(Pe) - 1)
    err = maximum(abs.(m.phi .- exact))
    # N=20 cannot fully resolve the outlet boundary layer (thickness ~1/Pe = 0.04 < h = 0.05);
    # the scheme is 2nd-order: max error 0.16 (N=20), 0.044 (N=40), 0.012 (N=80)
    check("1D central vs analytic (loose)", err < 0.2, @sprintf("max err = %.2e", err))
end

safely("Advection_Diffusion_OneDim_Upwind") do
    m = load("Advection_Diffusion_OneDim_Upwind.jl")
    ref = direct_1d(5, 2.5, 1.0, 1.0, 0.1, :upwind)
    err = maximum(abs.(m.phi .- ref))
    check("1D upwind vs direct solve", err < 1e-4, @sprintf("max err = %.2e", err))
    check("1D upwind bounded in [0,1]", all(0 .<= m.phi .<= 1), "")
end

# ---------- 1D conduction: exact solution is linear ----------

safely("Demo_1DProblem") do
    m = load("Demo_1DProblem.jl")
    err = maximum(abs.(m.T .- collect(m.x)))
    check("1D conduction vs linear profile", err < 1e-4, @sprintf("max err = %.2e", err))
end

# ---------- 2D convection-diffusion: direct solve of the same discrete system ----------
# Coefficients: W = D+F, S = D+Fy (boundaries 2D+F, 2D+Fy); E,N = D (boundary 2D); aP = sum of the four.

function direct_2d(N, L, rho, gamma, u, v, phiL, phiR, phiT, phiB)
    h = L / N; D = gamma / h; F = rho * u; Fy = rho * v
    idx(i, j) = (j - 1) * N + i
    A = zeros(N^2, N^2); b = zeros(N^2)
    for j in 1:N, i in 1:N
        k = idx(i, j); aP = 0.0
        # west
        if j > 1; A[k, idx(i,j-1)] = -(D+F); aP += D+F
        else b[k] += (2D+F)*phiL; aP += 2D+F end
        # east
        if j < N; A[k, idx(i,j+1)] = -D; aP += D
        else b[k] += 2D*phiR; aP += 2D end
        # north (row i-1)
        if i > 1; A[k, idx(i-1,j)] = -D; aP += D
        else b[k] += 2D*phiT; aP += 2D end
        # south (row i+1)
        if i < N; A[k, idx(i+1,j)] = -(D+Fy); aP += D+Fy
        else b[k] += (2D+Fy)*phiB; aP += 2D+Fy end
        A[k, k] = aP
    end
    return reshape(A \ b, N, N)
end

safely("TwoDim_Conv_Diff") do
    m = load("TwoDim_Conv_Diff.jl")
    ref = direct_2d(25, 1.0, 1.0, 0.0001, 5.0, 1.0, 1.0, 0.0, 1.0, 0.0)
    err = maximum(abs.(m.phi .- ref))
    check("2D conv-diff vs direct solve", err < 1e-3, @sprintf("max err = %.2e", err))
    check("2D conv-diff bounded in [0,1]", all(-1e-6 .<= m.phi .<= 1 + 1e-6), "")
end

# ---------- interactive / printing scripts (run as subprocesses with canned input) ----------

safely("DemoScript") do
    x0, h = 1.5, 0.1
    lines = run_script("DemoScript.jl"; input="$x0\n$h\n")
    fd = (f(x0 + h) - f(x0)) / h
    bd = (f(x0) - f(x0 - h)) / h
    cd = (f(x0 + h) - f(x0 - h)) / (2h)
    val(tag) = parse(Float64, lines[findfirst(==(tag), lines) + 1])
    ok = isapprox(val("Forward difference approximate"), fd; rtol=1e-8) &&
         isapprox(val("Backward difference approximate"), bd; rtol=1e-8) &&
         isapprox(val("Central difference approximate"), cd; rtol=1e-8)
    check("DemoScript finite differences", ok, "")
    # central difference of a cubic: f' + f'''h^2/6 with f''' = -24
    check("DemoScript central = f' - 4h^2", isapprox(cd, -12x0^2 + 14x0 - 3 - 4h^2; rtol=1e-8), "")
end

safely("Programming_Demo") do
    x0, h = 1.5, 0.1
    lines = run_script("Programming_Demo.jl"; input="$x0\n$h\n")
    ok = isapprox(parse_show(lines, "FDM"), (f(x0+h) - f(x0))/h; rtol=1e-8) &&
         isapprox(parse_show(lines, "BDM"), (f(x0) - f(x0-h))/h; rtol=1e-8) &&
         isapprox(parse_show(lines, "CDM"), (f(x0+h) - f(x0-h))/(2h); rtol=1e-8)
    check("Programming_Demo finite differences", ok, "")
end

safely("DemoScript_2") do
    lines = run_script("DemoScript_2.jl")
    expected = vcat(["Array Definition", "Array Calling", "-4.5", "4.5",
                     "for statement - default"], string.(1:10),
                    ["for statement - modified spacing"], string.(1:2:9),
                    ["if statement", "No", "while statement", "1", "2", "3"])
    check("DemoScript_2 exact output", lines == expected, "")
end

safely("Programming_Demo_2") do
    a = run_script("Programming_Demo_2.jl"; input="10\n")
    b = run_script("Programming_Demo_2.jl"; input="20\n")
    ok = count(==("For loop is executed"), a) == 5 &&
         "i is between 5 and 15" in a && !("i is greater than 2" in a) &&
         "i is not between 5 and 15" in b
    check("Programming_Demo_2 branches", ok, "")
end

# ---------- summary ----------

npass = count(last, results)
println("\n$npass / $(length(results)) checks passed")
exit(npass == length(results) ? 0 : 1)
