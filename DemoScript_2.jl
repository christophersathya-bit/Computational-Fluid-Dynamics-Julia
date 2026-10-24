# Equivalent of DemoScript_2.py

# Definition of arrays
println("Array Definition")
T1 = zeros(10)
T2 = ones(15)
T3 = Float64[1, 2, 3, 4]
T4 = collect(range(-4.5, 4.5, length=5))

# Indexing - Julia is 1-based
println("Array Calling")
println(T4[1])          # first element has index 1 in Julia
println(T4[end])        # last element

# for loop - default (Julia ranges include the end point)
println("for statement - default")
for i in 1:10
    println(i)
end

println("for statement - modified spacing")
for i in 1:2:9
    println(i)
end

# if statement
println("if statement")
if length(T4) == 4
    println("Yes")
else
    println("No")
end

# while statement
println("while statement")
counter = 0
while counter < 3
    global counter += 1
    println(counter)
end
