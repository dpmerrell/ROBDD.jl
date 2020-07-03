# BDD.jl
# 2020-07-02
# David Merrell


struct BDDNode end

struct DecisionNode <: BDDNode
    name::Symbol
    hi::BDDNode
    lo::BDDNode
end

struct TerminalNode <: BDDNode
    val::Bool
end



