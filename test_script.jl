
using ROBDD

my_bddtab = ROBDDTable([:x1,:x2,:x3])

myexpr = Meta.parse("!x1 | (x2 & !x2) | x3")

u = build_robdd(my_bddtab, myexpr)
println(my_bddtab)
println(u)

constraints = Dict(:x2 => true)
println("RESTRICT: ", constraints)
res_u = restrict(my_bddtab, u, constraints)
println(res_u)

println("CLEAN TABLE")
new_table, clean_idx = clean_table(my_bddtab, res_u)
println(new_table)
println(clean_idx)
