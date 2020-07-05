
using ROBDD

my_bddtab = ROBDDTable([:x1,:x2,:x3])

myexpr = Meta.parse("!x1 | (x2 & !x2) | x3")

u = build_robdd(my_bddtab, myexpr)
println(show(my_bddtab))
println(u)

constraints = Dict(:x2 => true)
println("RESTRICT: ", constraints)
res_u = restrict(my_bddtab, u, constraints)

println(res_u)
