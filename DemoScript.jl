# Equivalent of DemoScript.py

# Evaluate polynomial with coefficients in descending order (like np.polyval)
polyval(p, x) = foldl((acc, c) -> acc * x + c, p; init=zero(x))

function prompt(msg)
    print(msg)
    return parse(Float64, readline())
end

x0 = prompt("Enter the reference location: ")
h = prompt("Enter the spacing: ")

poly_p = [-4, 7, -3, 9]

F_x0 = polyval(poly_p, x0)
F_x0_plus_h = polyval(poly_p, x0 + h)
F_x0_minus_h = polyval(poly_p, x0 - h)

FDM = (F_x0_plus_h - F_x0) / h
println("Forward difference approximate")
println(FDM)

BDM = (F_x0 - F_x0_minus_h) / h
println("Backward difference approximate")
println(BDM)

CDM = (F_x0_plus_h - F_x0_minus_h) / (2h)
println("Central difference approximate")
println(CDM)
