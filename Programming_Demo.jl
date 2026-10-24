# Equivalent of Programming_Demo.m

polyval(p, x) = foldl((acc, c) -> acc * x + c, p; init=zero(x))

# Taking input from the user
print("Enter the value of x0: ")
x0 = parse(Float64, readline())
print("Enter the value of h: ")
h = parse(Float64, readline())

# Algebraic operations
a = x0 + h
b = x0 - h
c = x0 * h
d = x0 / h

# Polynomial evaluation
Fx = -4*x0^3 + 7*x0^2 - 3*x0 + 9
poly_p = [-4, 7, -3, 9]
F_X0 = polyval(poly_p, x0)
F_X0_Plus_H = polyval(poly_p, x0 + h)
F_X0_Minus_H = polyval(poly_p, x0 - h)

# Forward differencing
FDM = (F_X0_Plus_H - F_X0) / h
@show FDM

# Backward differencing
BDM = (F_X0 - F_X0_Minus_H) / h
@show BDM

# Central differencing
CDM = (F_X0_Plus_H - F_X0_Minus_H) / (2h)
@show CDM

# Displaying the output
println("Central differencing based estimate")
@show CDM
