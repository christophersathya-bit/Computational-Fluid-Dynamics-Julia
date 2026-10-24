# Equivalent of Programming_Demo_2.m

# Definition of arrays/vectors
T = [0, 0, 0, 0, 1]             # manual vector
T1 = zeros(3, 3)
T1_1 = zeros(3, 1)
T2 = ones(15, 1)
T2_2 = ones(1, 15)
T3 = [1, 2, 3, 4]
T4 = collect(range(-4.5, 4.5, length=5))

# Indexing (1-based, same as MATLAB)
T4[2] = 100                     # overwrite 2nd element

T_initial = zeros(1, 5)
T_initial[5] = 1

# for loop
for i in 1:5
    @show i
    println("For loop is executed")
end

# if statement
print("Enter a number: ")
i = parse(Float64, readline())
if 5 <= i <= 15
    println("i is between 5 and 15")
else
    println("i is not between 5 and 15")
end

# while statement (never executes, as in the original)
i = 1
while i > 2
    println("i is greater than 2")
    global i -= 1
    @show i
end

# The original MATLAB file contained pseudo-code here (not runnable):
#
# while numerical_error > required_epsilon
#     we want to keep doing the calculations (for)
#     repetitive iteration at different points
# end

# Data types
# Number - integer (Int), floating point (Float64)
# String - a single letter (Char), word/sentence (String)
