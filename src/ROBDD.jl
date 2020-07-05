# ROBDD.jl
# David Merrell
# (c) 2020-07-03
#
# A simple Julia module for building and working with ROBDDs.

module ROBDD

# The simple building block of our ROBDDs.
struct BDDNode
    var::Int64
    lo::Int64
    hi::Int64
end

# We define a `ROBDDContext` that 
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

struct ROBDDContext
    order::Vector{Symbol}
    order_lookup::Dict{Symbol,Int64}
    T::NodeTable
    H::LookupTable
end

ROBDDContext(order::Vector{Symbol}) = ROBDDContext(order, 
                                                   Dict(name=>idx for (idx,name) in enumerate(order)), 
                                                   NodeTable(length(order)), 
                                                   LookupTable())

# Some essential functionality for the ROBDDContext 
# and its members.

import Base: getindex, setindex!

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

function setindex!(H::LookupTable, idx::Int64, node::BDDNode)
    setindex!(H.d, idx, node)
end

function get_or_add(ctxt::ROBDDContext, node::BDDNode)
    if node.lo == node.hi
        return node.lo
    elseif member(ctxt.H, node)
        return ctxt.H[node]
    else
        u = add!(ctxt.T, node)
        ctxt.H[node] = u
        return u
    end
end

var(ctxt::ROBDDContext, idx::Int64) = ctxt.T[idx].var
lo(ctxt::ROBDDContext, idx::Int64) = ctxt.T[idx].lo
hi(ctxt::ROBDDContext, idx::Int64) = ctxt.T[idx].hi

"""
Apply a unary operation `op` to a BDD
"""
function apply(ctxt::ROBDDContext, op::Function, idx::Int64;
               memo::Dict{Tuple{Function,Vararg{Int64}},Int64}=Dict{Tuple{Function,Vararg{Int64}},Int64}())
    
    function rec_app(u1)
        if haskey(memo, (op,u1))
            return memo[(op,u1)]
        elseif in(u1, (0,1))
            u = Int64(op(Bool(u1)))
        else
            u = get_or_add(ctxt, 
                           BDDNode(var(ctxt,u1), 
                                   rec_app(lo(ctxt,u1)),
                                   rec_app(hi(ctxt,u1))
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
function apply(ctxt::ROBDDContext, op::Function, idx1::Int64, idx2::Int64; 
               memo::Dict{Tuple{Function,Vararg{Int64}},Int64}=Dict{Tuple{Function,Vararg{Int64}},Int64}())

    function rec_app(u1, u2)
        if haskey(memo, (op,u1,u2))
            return memo[(op,u1,u2)]
        elseif in(u1, (0,1)) & in(u2, (0,1))
            u = Int64(op(Bool(u1), Bool(u2)))
        elseif var(ctxt, u1) == var(ctxt, u2)
            u = get_or_add(ctxt, 
                           BDDNode(var(ctxt,u1), 
                                   rec_app(lo(ctxt,u1), lo(ctxt,u2)), 
                                   rec_app(hi(ctxt,u1), hi(ctxt,u2))
                                   )
                           )
        elseif var(ctxt, u1) < var(ctxt, u2)
            u = get_or_add(ctxt,
                           BDDNode(var(ctxt,u1),
                                   rec_app(lo(ctxt, u1), u2),
                                   rec_app(hi(ctxt, u1), u2)
                                   )
                           )
        else
            u = get_or_add(ctxt,
                           BDDNode(var(ctxt,u2),
                                   rec_app(u1, lo(ctxt, u2)),
                                   rec_app(u1, hi(ctxt, u2))
                                   )
                           )
        end
        memo[(op,u1,u2)] = u
        return u
    end
    
    return rec_app(idx1, idx2)
end

function build_robdd(ctxt::ROBDDContext, bool_expr::Union{Bool,Symbol,Expr})
    memo = Dict{Tuple{Function,Vararg{Int64}},Int64}()
    
    function rec_build(my_expr)
        # Terminals
        if typeof(my_expr) == Bool
            return Int64(my_expr)
        elseif typeof(my_expr) == Symbol
            u = get_or_add(ctxt, 
                           BDDNode(ctxt.order_lookup[my_expr],0,1)
                           )
        # Boolean op nodes
        else
            op = eval(my_expr.args[1])
            sub_exprs = my_expr.args[2:end]
            sub_bdds = Vector{Int64}(undef,length(sub_exprs))
            for (i,se) in enumerate(sub_exprs)
                sub_bdds[i] = rec_build(se)
            end
            u = apply(ctxt, op, sub_bdds...; memo=memo)
        end
        return u
    end
    
    return rec_build(bool_expr)
end


end # END MODULE

#my_ctxt = ROBDDContext([:x1,:x2])
#
#myexpr = Meta.parse("!x1 | x2")
#
#build_robdd(my_ctxt, myexpr)
#
#new_ctxt = ROBDDContext([:x2,:x1])
#
#build_robdd(new_ctxt, myexpr)
#
#new_ctxt.T[5]


