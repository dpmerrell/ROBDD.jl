# ROBDD.jl
# David Merrell
# (c) 2020-07-03
#
# A simple Julia module for building and working with ROBDDs.

module ROBDD


export ROBDDTable, build_robdd, apply, restrict, satcount, anysat, allsat, clean_table

# Some essential functionality for the ROBDDTable 
# and its members.
import Base: getindex, setindex!, show, iterate

using LRUCache

# The simple building block of our ROBDDs.
struct BDDNode
    var::UInt32
    lo::UInt32
    hi::UInt32
end


# We define a `ROBDDTable` that 
# (1) tracks the variables we're working with, and their ordering; and
# (2) maintains some hash tables for bookkeeping

mutable struct NodeTable
    d::Dict{UInt32,BDDNode}
    latest::UInt32
end


function NodeTable(n::UInt32)
    d = Dict{UInt32,BDDNode}()
    d[0] = BDDNode(n+1,0,0)
    d[1] = BDDNode(n+1,1,1)
    return NodeTable(d,1)
end


mutable struct LookupTable
    d::Dict{BDDNode,UInt32}
end


LookupTable() = LookupTable(Dict{BDDNode,UInt32}())


mutable struct ROBDDTable
    order::Vector{Symbol}
    order_lookup::Dict{Symbol,UInt32}
    T::NodeTable
    H::LookupTable
end


ROBDDTable(order::Vector{Symbol}) = ROBDDTable(order, 
                                               Dict(name=>idx for (idx,name) in enumerate(order)), 
                                               NodeTable(length(order)), 
                                               LookupTable()
                                              )



function show(io::IO, bddtab::ROBDDTable)
    out_str = string("(0) false\n(1) true")
    for i=2:bddtab.T.latest
        v = varname(bddtab,i)
        out_str = string(out_str, "\n","(",i,") ", 
                         varname(bddtab,i), " ", 
                         lo(bddtab,i), " ",
                         hi(bddtab,i))
    end
    print(io, out_str)

    return 
end


# Basic operations for the node table
function add!(T::NodeTable, node::BDDNode)
    T.latest = T.latest + 1
    T.d[T.latest] = node
    return T.latest
end

getindex(T::NodeTable, idx::UInt32) = T.d[idx]

# Basic operations for the "inverse" lookup table
member(H::LookupTable, node::BDDNode) = haskey(H.d, node)
getindex(H::LookupTable, node::BDDNode) = H.d[node]
getindex(bddtab::ROBDDTable, node::BDDNode) = bddtab.H[node]


function setindex!(H::LookupTable, idx::UInt32, node::BDDNode)
    setindex!(H.d, idx, node)
end


function get_or_add(bddtab::ROBDDTable, node::BDDNode)
    if node.lo == node.hi
        return node.lo
    elseif member(bddtab.H, node)
        return bddtab.H[node]
    else
        u = add!(bddtab.T, node)
        bddtab.H[node] = u
        return u
    end
end


getindex(bddtab::ROBDDTable, idx::UInt32) = bddtab.T[idx]
var(bddtab::ROBDDTable, idx::UInt32) = bddtab[idx].var
varname(bddtab::ROBDDTable, idx::UInt32) = bddtab.order[bddtab[idx].var]
lo(bddtab::ROBDDTable, idx::UInt32) = bddtab[idx].lo
hi(bddtab::ROBDDTable, idx::UInt32) = bddtab[idx].hi


"""
Apply a unary operation `op` to a BDD
"""
function apply(bddtab::ROBDDTable, op::Function, idx::UInt32;
               memo::AbstractDict{Tuple{Function,Vararg{UInt32}},UInt32}=LRU{Tuple{Function,Vararg{UInt32}},UInt32}(maxsize=1000000, by=sizeof)
              # memo::Dict{Tuple{Function,Vararg{UInt32}},UInt32}=Dict{Tuple{Function,Vararg{UInt32}},UInt32}())
              )
   
    function rec_app(u1)
        if haskey(memo, (op,u1))
            return memo[(op,u1)]
        elseif in(u1, (0,1))
            u = UInt32(op(Bool(u1)))
        else
            u = get_or_add(bddtab, 
                           BDDNode(var(bddtab,u1), 
                                   rec_app(lo(bddtab,u1)),
                                   rec_app(hi(bddtab,u1))
                                  )
                           )
        end
        memo[(op,u1)] = u
        return u
    end
        
    return rec_app(idx)
end


"""
Apply a binary operation `op` to two BDDs 
"""
function apply(bddtab::ROBDDTable, op::Function, idx1::UInt32, idx2::UInt32; 
               memo::AbstractDict{Tuple{Function,Vararg{UInt32}},UInt32}=LRU{Tuple{Function,Vararg{UInt32}},UInt32}(maxsize=1000000, by=sizeof)
              # memo::Dict{Tuple{Function,Vararg{UInt32}},UInt32}=Dict{Tuple{Function,Vararg{UInt32}},UInt32}()
              )


    function rec_app(u1, u2)
        if haskey(memo, (op,u1,u2))
            return memo[(op,u1,u2)]
        # Added some more cases to prune the recursion
        elseif (op == |) & (1 in (u1, u2))
            u = 1
        elseif (op == &) & (0 in (u1, u2))
            u = 0
        elseif in(u1, (0,1)) & in(u2, (0,1))
            u = UInt32(op(Bool(u1), Bool(u2)))
        elseif var(bddtab, u1) == var(bddtab, u2)
            u = get_or_add(bddtab, 
                           BDDNode(var(bddtab,u1), 
                                   rec_app(lo(bddtab,u1), lo(bddtab,u2)), 
                                   rec_app(hi(bddtab,u1), hi(bddtab,u2))
                                   )
                           )
        elseif var(bddtab, u1) < var(bddtab, u2)
            u = get_or_add(bddtab,
                           BDDNode(var(bddtab,u1),
                                   rec_app(lo(bddtab, u1), u2),
                                   rec_app(hi(bddtab, u1), u2)
                                   )
                           )
        else
            u = get_or_add(bddtab,
                           BDDNode(var(bddtab,u2),
                                   rec_app(u1, lo(bddtab, u2)),
                                   rec_app(u1, hi(bddtab, u2))
                                   )
                           )
        end
        memo[(op,u1,u2)] = u
        return u
    end
    
    return rec_app(idx1, idx2)
end


"""
Build an ROBDD from a boolean expression `bool_expr`,
within the context defined by ROBDDTable `bddtab`.

Optional `memo` kwarg allows a persistent memoization cache 
to be used by multiple calls to `build_robdd`.
"""
function build_robdd(bddtab::ROBDDTable, bool_expr::Union{Bool,Symbol,Expr};
                     memo::AbstractDict{Tuple{Function,Vararg{UInt32}},UInt32}=LRU{Tuple{Function,Vararg{UInt32}},UInt32}(maxsize=1000000, by=sizeof)
                    )
    
    function rec_build(my_expr)
        # Terminals
        if typeof(my_expr) == Bool
            return UInt32(my_expr)
        elseif typeof(my_expr) == Symbol
            u = get_or_add(bddtab, 
                           BDDNode(bddtab.order_lookup[my_expr],0,1)
                           )
        # Boolean op nodes
        else
            u = -1
            op = eval(my_expr.args[1])
            sub_exprs = my_expr.args[2:end]
            sub_bdds = Vector{UInt32}(undef,length(sub_exprs))
            for (i,se) in enumerate(sub_exprs)
                sub_bdds[i] = rec_build(se)
                # Add some checks for lazy evaluation
                if (op == |) & (sub_bdds[i] == 1)
                    u = 1
                    break
                elseif (op == &) & (sub_bdds[i] == 0)
                    u = 0
                    break
                end
            end
            if u == -1
                u = apply(bddtab, op, sub_bdds...; memo=memo)
            end
        end
        return u
    end
    
    return rec_build(bool_expr)
end


"""
Restrict a ROBDD via `assignments`,
a dictionary of variable assignments.

Optional `memo` kwarg allows a persistent memoization cache 
to be used by multiple calls to `restrict`.
"""
function restrict(bddtab::ROBDDTable, bdd_idx::UInt32, assignments::Dict{Symbol,Bool};
                  memo::AbstractDict{Tuple{Function,Vararg{UInt32}},UInt32}=LRU{Tuple{Function,Vararg{UInt32}},UInt32}(maxsize=1000000, by=sizeof)
                 # memo::Dict{Tuple{UInt32,UInt32,Bool},UInt32}=Dict{Tuple{UInt32,UInt32,Bool},UInt32}()
                  )
    
    function rec_res(u, j, b)
        if haskey(memo, (u,j,b))
            return memo[(u,j,b)]
        elseif var(bddtab, u) > j
            return u
        elseif var(bddtab, u) < j
            u_new = get_or_add(bddtab, BDDNode(var(bddtab, u),
                                               rec_res(lo(bddtab, u), j, b),
                                               rec_res(hi(bddtab, u), j, b)
                                               )
                               )
        elseif !b
            u_new = rec_res(lo(bddtab,u), j, b)
        else
            u_new = rec_res(hi(bddtab,u), j, b)
        end
        memo[(u,j,b)] = u_new
        return u_new
    end

    for (j, b) in pairs(assignments)
        bdd_idx = rec_res(bdd_idx, bddtab.order_lookup[j], b)
    end

    return bdd_idx
end


"""
Count the satisfying assignments of an ROBDD.
"""
function satcount(bddtab::ROBDDTable, idx::UInt32;
                  memo::AbstractDict{Tuple{Function,Vararg{UInt32}},UInt32}=LRU{Tuple{Function,Vararg{UInt32}},UInt32}(maxsize=1000000, by=sizeof)
                  #memo::Dict{UInt32,UInt32}=Dict{UInt32,UInt32}())
                  )

    function rec_count(u)
        if haskey(memo, u)
            return memo[u]
        elseif u == 0
            count = 0
        elseif u == 1
            count = 1
        else
            l = lo(bddtab,u)
            h = hi(bddtab,u)
            lcount = rec_count(l)
            hcount = rec_count(h)
            res = lcount * 2^(var(bddtab,l) - var(bddtab,u) - 1) + hcount * 2^(var(bddtab,r) - var(bddtab,u) - 1)
        end
        return res
    end

    return rec_count(idx)
end


####################################################
# Helper functions for SAT iteration
####################################################
function state_to_assignment(bddtable::ROBDDTable, 
                             state::Vector{Pair{UInt32,Bool}})
    assignment = [bddtable.order[var(bddtable,n)] => b for (n,b) in state]
    assignment = Dict{Symbol,Bool}(assignment)
    return assignment
end


function dfs_sat_first(bddtab::ROBDDTable, idx::UInt32)
    if idx == 0
        error("unsat")
    elseif idx == 1
        result = Vector{Pair{UInt32,Bool}}()
    elseif lo(bddtab, idx) == 0
        result = dfs_sat_first(bddtab, hi(bddtab,idx))
        pushfirst!(result, idx => true)
    else
        result = dfs_sat_first(bddtab, lo(bddtab,idx))
        pushfirst!(result, idx => false)
    end
    return result
end


function is_sat(bddtab::ROBDDTable, state)
    (idx,b) = state[end]
    if b & hi(bddtab, idx) == 1
        return true
    elseif !b & lo(bddtab, idx) == 1
        return true
    else
        return false
    end
end


function increment_ancestor(bddtab::ROBDDTable, state)
    n = length(state)
    for i = n:-1:1
        (idx,b) = pop!(state)
        if !b
            push!(state, idx => true)
            break
        end
    end
    return state
end


function dfs_next_sat(bddtab::ROBDDTable, state)

    while true
        (idx,b) = state[end]
        if !b & !in(lo(bddtab,idx), (0,1))
            # lo child not terminal --> add lo child
            push!(state, lo(bddtab,idx) => false)
        elseif b & !in(hi(bddtab,idx), (0,1))
            # hi child not terminal --> add hi child
            push!(state, hi(bddtab,idx) => false)
        else
            # change self or nearest ancestor from false -> true
            state = increment_ancestor(bddtab, state)
        end

        if length(state) == 0 
            break
        elseif is_sat(bddtab, state)
            break
        end
    end
    return state
end


"""
Return an arbitrary satisfying assignment for an ROBDD.
The method is deterministic -- it will return the same value
when called multiple times.
"""
function anysat(bddtab::ROBDDTable, idx::UInt32)
    assignments = dfs_sat_first(bddtab, idx)
    assignments = [bddtab.order[var(bddtab,idx)] => b for (idx,b) in assignments]
    return Dict{Symbol,Bool}(assignments)
end


##########################################
# AllSat iterator
##########################################

"""
An iterator over the set of satisfying assignments.
Performs DFS to find all paths from the root node 
to the (1) node.

Note that any unassigned variables may take arbitrary values.
Hence, an assignment with n unassigned variables actually 
represents 2^n assignments.
"""
struct AllSat
    bddtable::ROBDDTable
    idx::UInt32
end

allsat(bddtab::ROBDDTable, idx::UInt32) = AllSat(bddtab,idx)


function iterate(iter::AllSat)
    try 
        state = dfs_sat_first(iter.bddtable, iter.idx)
        soln = state_to_assignment(iter.bddtable, state)
        return (soln, state)
    catch 
        return nothing
    end
end


function iterate(iter::AllSat, state::Vector{Pair{UInt32,Bool}})
    state = dfs_next_sat(iter.bddtable, state)
    if length(state) == 0
        return nothing
    else
        return (state_to_assignment(iter.bddtable, state), state)
    end
end


"""
Create a new ROBDDTable encoding the same ROBDD
as `bddtab[idx]`, but with all of the "unused" nodes
removed.

Returns a tuple: (new_table::ROBDDTable, new_idx::UInt32)
"""
function clean_table(bddtab::ROBDDTable, idx::UInt32)
    
    new_table = ROBDDTable(bddtab.order)
    old2new = Dict{UInt32,UInt32}(0=>0,1=>1)

    function rec_clean(old_u)
        if haskey(old2new, old_u)
            return old2new[old_u]
        else
            new_lo = rec_clean(lo(bddtab, old_u))
            new_hi = rec_clean(hi(bddtab, old_u))
            new_u = get_or_add(new_table, 
                               BDDNode(var(bddtab, old_u),
                                       new_lo, new_hi
                                      )
                               )
            old2new[old_u] = new_u
            return new_u
        end
    end

    new_idx = rec_clean(idx)
    return new_table, new_idx
end


end # END MODULE


