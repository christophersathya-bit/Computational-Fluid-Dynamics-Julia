# Computational Fluid Dynamics in Julia – a beginner's path

Learn the core ideas of **Computational Fluid Dynamics (CFD)** by running small, readable solvers written in **Julia**. You start with "how does a computer approximate a derivative?" and end with a 2D convection–diffusion solver. Every result is checked against an exact or reference answer.

> **Who is this for?** Students and engineers who know basic calculus and have written a little code in any language. No previous CFD or Julia experience is needed.

---

## Contents

1. [What you will learn](#1-what-you-will-learn)
2. [Install Julia (step by step)](#2-install-julia-step-by-step)
3. [Get the code and run your first script](#3-get-the-code-and-run-your-first-script)
4. [Learning path – run the scripts in this order](#4-learning-path--run-the-scripts-in-this-order)
5. [The physics and the numerics, explained](#5-the-physics-and-the-numerics-explained)
6. [Validation – how we know the answers are right](#6-validation--how-we-know-the-answers-are-right)
7. [Exercises to try](#7-exercises-to-try)
8. [Troubleshooting](#8-troubleshooting)
9. [Glossary](#9-glossary)
10. [Further reading](#10-further-reading)

---

## 1. What you will learn

| Concept | Where you meet it |
| --- | --- |
| Forward, backward and central **finite differences** | `Programming_Demo.jl`, `DemoScript.jl` |
| Julia arrays, loops and if/while statements | `Programming_Demo_2.jl`, `DemoScript_2.jl` |
| Turning a differential equation into a **system of algebraic equations** | `Demo_1DProblem.jl` |
| **Iterative solvers** (Jacobi) and **convergence** (residuals) | all solver scripts |
| The **finite volume method**: cells, faces, fluxes, boundary conditions | `Advection_Diffusion_OneDim*.jl` |
| The **Péclet number**, **central vs upwind** schemes and **stability** | `Advection_Diffusion_OneDim*.jl` |
| Extending a solver to **2D**, including edge and corner boundaries | `TwoDim_Conv_Diff.jl` |
| **Verification**: comparing with analytic and direct solutions, grid refinement | `validate.jl`, `reference_check.py` |

---

## 2. Install Julia (step by step)

1. Go to **https://julialang.org/downloads/** and install Julia **1.9 or newer**.
   - **Windows:** the easiest option is the Microsoft Store app "Julia" (it installs `juliaup`). You can also run `winget install julia -s msstore` in a terminal.
   - **macOS / Linux:** run `curl -fsSL https://install.julialang.org | sh` in a terminal.
2. Open a new terminal and type `julia`. You should see the Julia prompt:
   ```
   julia>
   ```
3. Install the plotting package (this needs to be done only once):
   ```julia
   julia> using Pkg
   julia> Pkg.add("Plots")
   ```
   The first install and the first plot can take a few minutes, because Julia compiles packages the first time. Later runs are much faster.

**Recommended editor:** [VS Code](https://code.visualstudio.com/) with the **Julia extension**. It lets you run a file with one click and shows plots in a side panel.

---

## 3. Get the code and run your first script

**Download the code**

- Click the green **Code** button on this page, then **Download ZIP**, and unzip it, **or**
- if you use git: `git clone https://github.com/christophersathya-bit/Computational-Fluid-Dynamics-Julia.git`

**Run a script (the beginner-friendly way)**

1. Open a terminal **inside the downloaded folder**.
2. Start Julia with this project's settings:
   ```
   julia --project=.
   ```
3. At the `julia>` prompt, run a script with `include`:
   ```julia
   julia> include("Demo_1DProblem.jl")
   ```
   The script prints the number of iterations and opens a plot window.

> **Why `include` instead of `julia script.jl`?** When you run a file from the terminal, Julia exits as soon as the script ends, so plot windows close immediately. Running it with `include` inside the Julia prompt keeps the plots on screen.

**Scripts that ask for input**

`Programming_Demo.jl`, `DemoScript.jl` and `Programming_Demo_2.jl` ask you to type numbers. For example, enter `1.5` for x₀ and `0.1` for h, pressing Enter after each one.

**Run everything at once (Linux / macOS / Git Bash)**

```bash
./run_all.sh
```

This installs packages, runs the validation suite and runs every script without opening windows. Output goes to the `logs/` folder.

---

## 4. Learning path – run the scripts in this order

| Step | Script | What to look at |
| --- | --- | --- |
| 1 | `Programming_Demo.jl` / `DemoScript.jl` | Compare the forward, backward and central estimates of f′(x) with the exact value. Which one is closest? |
| 2 | `Programming_Demo_2.jl` / `DemoScript_2.jl` | Julia basics: arrays start at index **1**, `range`, loops, `if`, `while` |
| 3 | `Demo_1DProblem.jl` | Steady heat conduction in a rod. Watch the residual fall on the convergence plot. The answer should be a straight line. |
| 4 | `Advection_Diffusion_OneDim.jl` | Flow carries φ to the right, and diffusion spreads it. Note the steep drop near x = 1 (a **boundary layer**). |
| 5 | `Advection_Diffusion_OneDim_Upwind.jl` | Same problem on only 5 cells. The upwind scheme stays stable where the central scheme would fail. |
| 6 | `TwoDim_Conv_Diff.jl` | The same ideas in 2D: flow at an angle, with edges and corners to handle |
| 7 | `validate.jl` | Automatic checks of every script against reference answers |

---

## 5. The physics and the numerics, explained

### 5.1 Finite differences – approximating a derivative

A computer cannot take a limit, so we replace the derivative with a difference over a small step $h$:

| Name | Formula | Error |
| --- | --- | --- |
| Forward | $f'(x) \approx \dfrac{f(x+h) - f(x)}{h}$ | proportional to $h$ (first order) |
| Backward | $f'(x) \approx \dfrac{f(x) - f(x-h)}{h}$ | proportional to $h$ (first order) |
| Central | $f'(x) \approx \dfrac{f(x+h) - f(x-h)}{2h}$ | proportional to $h^2$ (second order) |

**Second order** means that halving $h$ makes the error about **4× smaller**. For the cubic in the demo ($f = -4x^3 + 7x^2 - 3x + 9$, $x_0 = 1.5$, $h = 0.1$) the exact derivative is −9.00. Forward gives −10.14, backward −7.94 and central −9.04.

### 5.2 Steady 1D heat conduction (`Demo_1DProblem.jl`)

The equation is $\dfrac{d^2T}{dx^2} = 0$, with $T(0) = 0$ and $T(1) = 1$. Replacing the second derivative with a central difference gives, at every interior point:

$$T_i = \tfrac{1}{2}\,(T_{i-1} + T_{i+1})$$

Each point is the average of its neighbours. The **Jacobi method** applies this update to every point again and again until the values stop changing. The **residual** (the total change in one sweep) measures how far we are from convergence. The script stops when it drops below $10^{-7}$.

### 5.3 1D advection–diffusion (`Advection_Diffusion_OneDim*.jl`)

$$\frac{d}{dx}(\rho u \phi) = \frac{d}{dx}\!\left(\Gamma \frac{d\phi}{dx}\right), \qquad \phi(0) = 1,\ \phi(L) = 0$$

- **Advection** ($\rho u \phi$): the flow carries φ along, here to the right with $u = 2.5$.
- **Diffusion** ($\Gamma\, d\phi/dx$): φ spreads out, with diffusion coefficient $\Gamma = 0.1$.

**Finite volume method.** The domain is split into $N$ cells. For each cell we balance what flows **in** and **out** through its two faces. This gives one equation per cell:

$$a_P\,\phi_P = a_W\,\phi_W + a_E\,\phi_E$$

with $F = \rho u$ (advection strength) and $D = \Gamma / h$ (diffusion strength):

| Scheme | $a_W$ | $a_E$ |
| --- | --- | --- |
| Central | $D + F/2$ | $D - F/2$ |
| Upwind | $D + F$ | $D$ |

At the boundaries the distance from the cell centre to the wall is only $h/2$, so the diffusion coefficient there doubles to $2D$.

**The Péclet number** compares advection with diffusion: $Pe = F/D = \rho u h / \Gamma$.

- **Central** differencing is more accurate, but only works if $Pe < 2$. Above that, $a_E$ becomes negative, the solution oscillates and the iteration can blow up.
- **Upwind** differencing takes information only from the direction the flow comes from. It is always stable and never oscillates, but it smears sharp changes ("numerical diffusion").

In this repo: 20 cells give $Pe = 1.25$, so central works. 5 cells give $Pe = 5$, so upwind is needed.

**Exact solution** (for comparison):

$$\phi(x) = 1 - \frac{e^{Pe_L\, x/L} - 1}{e^{Pe_L} - 1}, \qquad Pe_L = \frac{\rho u L}{\Gamma} = 25$$

### 5.4 2D convection–diffusion (`TwoDim_Conv_Diff.jl`)

The same balance in a unit square, with flow $u = 5$ (to the right) and $v = 1$ (upwards), and very weak diffusion ($\Gamma = 10^{-4}$). Each cell now has four neighbours (W, E, S, N):

$$a_P\,\phi_P = a_W\,\phi_W + a_E\,\phi_E + a_S\,\phi_S + a_N\,\phi_N + S_u$$

- Upwind in both directions: $a_W = D + F_x$ and $a_S = D + F_y$, because the flow comes **from** the west and the south.
- **Boundary values:** φ = 1 on the left and top walls, φ = 0 on the right and bottom walls. Cells next to a wall get an extra source term $S_u$.
- The code treats **9 cell types** separately: interior, 4 edges and 4 corners. This is the trickiest part to get right, so read it slowly.

---

## 6. Validation – how we know the answers are right

A CFD code can produce a nice-looking picture that is completely wrong. `validate.jl` checks every script against an answer that is known independently:

- **Iterative vs direct solve:** the same equations are also solved exactly with a matrix solve (`A \ b`). The two answers should agree to about $10^{-7}$. This checks the **iteration**.
- **Numerical vs analytic solution:** this checks the **discretisation**, i.e. how well the grid represents the real equation.
- **Boundedness:** φ must stay between 0 and 1 (no unphysical overshoots).
- **Grid refinement:** halving the cell size should make the error ~4× smaller for a second-order scheme.

```bash
julia --project=. validate.jl     # prints PASS/FAIL for each check, exit code 1 if any fail
```

### Results

![Validation results](validation_results.png)

| Check | Result |
| --- | --- |
| 1D central (N=20), iterative vs direct solve | max error 9.0e-08 |
| 1D upwind (N=5), iterative vs direct solve | max error 5.1e-08 |
| 1D conduction vs exact linear profile | max error 6.0e-07 |
| 2D convection–diffusion (25×25), iterative vs direct solve | max error 1.1e-07 |
| Upwind and 2D solutions bounded in [0, 1] | pass |
| Central difference of a cubic, error = −4h² | exact to 2.6e-14 |
| **Grid refinement of the central scheme vs analytic** | **0.160 → 0.044 → 0.012 → 0.003 (N = 20 → 160), second-order** |

**How to read these results**

- On the 20-cell grid the central scheme differs from the analytic solution by 0.16. The boundary layer near x = 1 (thickness ≈ 1/Pe = 0.04) is thinner than one cell (h = 0.05). The error falls about 4× each time the grid is doubled, which confirms second-order accuracy.
- On the 5-cell grid ($Pe = 5$) the central scheme is unstable and the iteration diverges. That is why the coarse-grid case uses **upwind**: it stays bounded, but adds numerical diffusion (error 0.20 vs analytic).
- In 2D, the low-φ region from the bottom wall (φ = 0) is carried up and to the right along the flow direction (slope v/u = 0.2). The inflow value φ = 1 from the left wall fills the rest of the domain. Diffusion is weak (Γ = 10⁻⁴), so the mixing layer between them stays thin.

The figure and table are produced by `reference_check.py`, an independent NumPy implementation of the same discrete equations. It cross-checks the Julia solvers without needing Julia (`python reference_check.py`, needs `numpy` and `matplotlib`).

---

## 7. Exercises to try

1. **Order of accuracy:** in `Programming_Demo.jl`, try $h$ = 0.1, 0.05 and 0.025. How fast does each error shrink?
2. **Break the central scheme:** in `Advection_Diffusion_OneDim.jl`, set `N = 5`. What happens to the residual, and why? (Hint: compute $Pe$.)
3. **Grid refinement:** run the central scheme with N = 10, 20, 40, 80 and plot the maximum error against N on log–log axes. Is the slope −2?
4. **Faster convergence:** replace Jacobi with **Gauss–Seidel** (use the new value of a neighbour as soon as it is computed). How many iterations do you save?
5. **Change the flow angle:** in `TwoDim_Conv_Diff.jl`, try `v = 5.0` (45°) and `v = 0.0` (horizontal flow) and compare the contour plots.
6. **Stronger diffusion:** set `gamma = 0.01` in 2D. Why does the sharp layer become wider?

---

## 8. Troubleshooting

| Problem | Fix |
| --- | --- |
| `ERROR: ArgumentError: Package Plots not found` | Run `using Pkg; Pkg.add("Plots")`, or start Julia with `julia --project=.` and run `Pkg.instantiate()` |
| The first run is very slow | Normal. Julia compiles packages on first use. The second run is fast. |
| The plot window flashes and closes | Run the script with `include("file.jl")` inside the Julia prompt (see section 3) |
| The script waits and does nothing | It is waiting for keyboard input. Type a number and press Enter. |
| Residual grows or prints `NaN` | The scheme is unstable. Check the Péclet number (central needs $Pe < 2$) or switch to upwind. |
| `./run_all.sh: Permission denied` | Run `chmod +x run_all.sh` once, or start it with `bash run_all.sh` |

---

## 9. Glossary

| Term | Meaning |
| --- | --- |
| **Discretisation** | Replacing a continuous equation with equations at a finite number of points or cells |
| **Grid / mesh** | The set of points or cells the domain is divided into. $h$ is the spacing between them. |
| **Finite difference** | Approximating derivatives with differences between neighbouring point values |
| **Finite volume** | Balancing fluxes in and out of each cell. It conserves quantities exactly, which is why most CFD codes use it. |
| **Advection / convection** | Transport of a quantity by the flow |
| **Diffusion** | Spreading of a quantity from high to low concentration |
| **Péclet number** | Ratio of advection to diffusion strength |
| **Jacobi iteration** | Repeatedly updating every unknown from its neighbours' old values until nothing changes |
| **Residual** | How much the solution changed in one iteration, used to judge convergence |
| **Boundedness** | The solution stays within physically sensible limits (here between 0 and 1) |
| **Numerical diffusion** | Artificial smearing caused by a low-order scheme such as upwind |
| **Order of accuracy** | How fast the error falls as $h$ shrinks: error ∝ $h^p$ for a $p$-th order scheme |

---

## 10. Further reading

- H. K. Versteeg & W. Malalasekera, *An Introduction to Computational Fluid Dynamics: The Finite Volume Method* (2nd ed.). The coefficient notation ($a_W$, $a_E$, $a_P$) used here follows this book.
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, a classic on finite volume methods.
- J. D. Anderson, *Computational Fluid Dynamics: The Basics with Applications*, a gentle introduction.
- The official Julia manual: https://docs.julialang.org

## Roadmap

- [ ] Lid-driven cavity (SIMPLE algorithm, staggered grid) validated against Ghia et al. (1982), Re = 100
- [ ] Gauss–Seidel / SOR and a sparse direct solver for faster convergence
- [ ] Higher-order schemes (QUICK, TVD limiters)
- [ ] Transient problems with explicit and implicit time stepping

## Repository structure

```
├── Programming_Demo.jl, DemoScript.jl          # finite-difference basics (interactive)
├── Programming_Demo_2.jl, DemoScript_2.jl      # Julia language basics
├── Demo_1DProblem.jl                           # 1D steady heat conduction
├── Advection_Diffusion_OneDim.jl               # 1D advection–diffusion, central scheme
├── Advection_Diffusion_OneDim_Upwind.jl        # 1D advection–diffusion, upwind scheme
├── TwoDim_Conv_Diff.jl                         # 2D convection–diffusion, upwind
├── validate.jl                                 # automated checks (PASS/FAIL)
├── reference_check.py                          # NumPy cross-check, produces validation_results.png
├── validation_results.png                      # validation figure
├── run_all.sh                                  # install + validate + run everything
└── Project.toml                                # Julia package list
```

---

*The model problems follow standard CFD coursework on finite-difference and finite-volume methods. This repository implements them in Julia, with an automated validation suite and a beginner's guide.*

**Christopher Sathiya Satchi**, Junior Researcher at CIIAE, Spain · [LinkedIn](https://www.linkedin.com/in/christophers828269107/) · [Google Scholar](https://scholar.google.com/citations?user=s9fiBJoAAAAJ)

Questions or suggestions? Open an **Issue** on this repository.
