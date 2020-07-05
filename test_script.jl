
using ROBDD

my_bddtab = ROBDDTable([:x1,:x2])

myexpr = Meta.parse("(x1 | !x1) | (x2 & !x2)")

u = build_robdd(my_bddtab, myexpr)
println(show(my_bddtab))
println(u)
