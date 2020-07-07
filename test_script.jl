# Load the module
using ROBDD

# Define your boolean variables
# (must be a vector of symbols)
ordered_variables = [:x1, :x2, :x3]

# Initiate a ROBDDTable. This datastructure
# encodes context information (the set of variables
# and their order); and *all* of the ROBDDs you 
# construct over them. Each entry is a node
# in a ROBDD.
my_bdd_table = ROBDDTable(ordered_variables)

# At first it only has two entries. 
# one for "false" and one for "true".
println(my_bdd_table)

# output:
#
# (0) false
# (1) true

# Construct the Boolean expression you're modeling
#my_expr = :(x1 | (x2 & !x2) | x3)
my_expr = :(x1 | x2 | x3)

# Build a ROBDD from the expression.
# The `build_robdd` function returns an 
# integer index. The index points to a particular
# entry of the ROBDDTable.
idx = build_robdd(my_bdd_table, my_expr)

# During construction, several new
# entries were added to your table.
# `build_robdd` tells you which entry
# contains the root node of your BDD.
println(my_bdd_table)
println("INDEX: ", idx)

# output:
#
# (0) false
# (1) true
# (2) x1 0 1
# (3) x2 0 1
# (4) x2 1 0
# (5) x3 0 1
# (6) x1 5 1
# INDEX: 6


# We provide an interface for
# (1) getting an arbitrary SAT solution and
# (2) enumerating all SAT solutions
#     (modulo the unassigned variables--these
#      may take any value.)
println("ANYSAT: ", anysat(my_bdd_table, idx))
println("ALLSAT:")
for sat in allsat(my_bdd_table, idx)
    println("\t",sat)
end

# output:
#
# ANYSAT: Dict{Symbol,Bool}(:x2 => 0,:x3 => 1,:x1 => 0)
# ALLSAT:
# 	Dict{Symbol,Bool}(:x2 => 0,:x3 => 1,:x1 => 0)
# 	Dict{Symbol,Bool}(:x2 => 1,:x1 => 0)
# 	Dict{Symbol,Bool}(:x1 => 1)


# We can "restrict" a ROBDD by assigning
# concrete values to a subset of its variables.
constraints = Dict(:x3 => false)
restricted_idx = restrict(my_bdd_table, idx, constraints)
println("RESTRICTIONS: ", constraints)
println(my_bdd_table)
println("RESTRICTED IDX: ", restricted_idx)

# output:
#
# RESTRICTIONS: Dict{Symbol,Bool}(:x3 => 0)
# (0) false
# (1) true
# (2) x1 0 1
# (3) x2 0 1
# (4) x2 1 0
# (5) x3 0 1
# (6) x1 5 1
# RESTRICTED IDX: 2

# During construction the ROBDDTable accumulates
# many entries. However, some of them are not actually used by the
# ROBDD of interest. We provide a function for "cleaning"
# the table, leaving only the descendants of a specified BDD.
new_table, clean_idx = clean_table(my_bdd_table, restricted_idx)
println("CLEANED TABLE:")
println(new_table)
println("CLEAN IDX: ", clean_idx)

# output:
#
# CLEANED TABLE:
# (0) false
# (1) true
# (2) x1 0 1
# CLEAN IDX: 2

# You should only clean a ROBDDTable when you don't intend to 
# construct any more ROBDDs. Otherwise, those extra ROBDDs
# may be useful for constructing other ROBDDs.

# ***MORE FEATURES TO COME***
