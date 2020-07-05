# ROBDD.jl
# David Merrell
# (c) 2020-07-03
#
# A simple Julia module for building and working with ROBDDs.

module ROBDD


export ROBDDTable, build_robdd, apply, restrict, clean_table 

# Some essential functionality for the ROBDDTable 
# and its members.
import Base: getindex, setindex!, show


# The simple building block of our ROBDDs.
struct BDDNode
    var::Int64
    lo::Int64
    hi::Int64
end


# We define a `ROBDDTable` that 
# (1) tracks the variables we're working with, and their ordering; and
# (2) maintains some hash tables for bookkeeping

mutable struct NodeTable
    d::Dict{Int64,BDDNode}
    latest::Int64
end


function NodeTable(n::Int64)
    d = Dict{Int64,BDDNode}()
    d[0] = BDDNode(n+1,0,0)
    d[1] = BDDNode(n+1,1,1)
    return NodeTable(d,1)
end


mutable struct LookupTable
    d::Dict{BDDNode,Int64}
end


LookupTable() = LookupTable(Dict{BDDNode,Int64}())


mutable struct ROBDDTable
    order::Vector{Symbol}
    order_lookup::Dict{Symbol,Int64}
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

getindex(T::NodeTable, idx::Int64) = T.d[idx]

# Basic operations for the "inverse" lookup table
member(H::LookupTable, node::BDDNode) = haskey(H.d, node)
getindex(H::LookupTable, node::BDDNode) = H.d[node]
getindex(bddtab::ROBDDTable, node::BDDNode) = bddtab.H[node]


function setindex!(H::LookupTable, idx::Int64, node::BDDNode)
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


getindex(bddtab::ROBDDTable, idx::Int64) = bddtab.T[idx]
var(bddtab::ROBDDTable, idx::Int64) = bddtab[idx].var
varname(bddtab::ROBDDTable, idx::Int64) = bddtab.order[bddtab[idx].var]
lo(bddtab::ROBDDTable, idx::Int64) = bddtab[idx].lo
hi(bddtab::ROBDDTable, idx::Int64) = bddtab[idx].hi


"""
Apply a unary operation `op` to a BDD
"""
function apply(bddtab::ROBDDTable, op::Function, idx::Int64;
               memo::Dict{Tuple{Function,Vararg{Int64}},Int64}=Dict{Tuple{Function,Vararg{Int64}},Int64}())
   
    function rec_app(u1)
        if haskey(memo, (op,u1))
            return memo[(op,u1)]
        elseif in(u1, (0,1))
            u = Int64(op(Bool(u1)))
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
function apply(bddtab::ROBDDTable, op::Function, idx1::Int64, idx2::Int64; 
               memo::Dict{Tuple{Function,Vararg{Int64}},Int64}=Dict{Tuple{Function,Vararg{Int64}},Int64}())


    function rec_app(u1, u2)
        if haskey(memo, (op,u1,u2))
            return memo[(op,u1,u2)]
        # Added some more cases to prune the recursion
        elseif (op == |) & (1 in (u1, u2))
            u = 1
        elseif (op == &) & (0 in (u1, u2))
            u = 0
        elseif in(u1, (0,1)) & in(u2, (0,1))
            u = Int64(op(Bool(u1), Bool(u2)))
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
                     memo::Dict{Tuple{Function,Vararg{Int64}},Int64}=Dict{Tuple{Function,Vararg{Int64}},Int64}()
                    )
    
    function rec_build(my_expr)
        # Terminals
        if typeof(my_expr) == Bool
            return Int64(my_expr)
        elseif typeof(my_expr) == Symbol
            u = get_or_add(bddtab, 
                           BDDNode(bddtab.order_lookup[my_expr],0,1)
                           )
        # Boolean op nodes
        else
            u = -1
            op = eval(my_expr.args[1])
            sub_exprs = my_expr.args[2:end]
            sub_bdds = Vector{Int64}(undef,length(sub_exprs))
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
function restrict(bddtab::ROBDDTable, bdd_idx::Int64, assignments::Dict{Symbol,Bool};
                  memo::Dict{Tuple{Int64,Int64,Bool},Int64}=Dict{Tuple{Int64,Int64,Bool},Int64}()
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
Create a new ROBDDTable encoding the same ROBDD
as `bddtab[idx]`, but with all of the "unused" nodes
removed.

Returns a tuple: (new_table::ROBDDTable, new_idx::Int64)
"""
function clean_table(bddtab::ROBDDTable, idx::Int64)
    
    new_table = ROBDDTable(bddtab.order)
    old2new = Dict{Int64,Int64}(0=>0,1=>1)

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


