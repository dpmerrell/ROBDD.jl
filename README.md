# BDD.jl

A Julia implementation of [Binary Decision Diagrams](https://en.wikipedia.org/wiki/Binary_decision_diagram).

> "...one of the only really fundamental data structures that came out in the last twenty-five years"
>
> -- Donald Knuth. [_Fun with Binary Decision Diagrams_](https://www.youtube.com/watch?v=SQE21efsf7Y). 2008


Some basic requirements:

* User should be able to manually construct a ROBDD by performing boolean operations on sub-ROBDDs
    - 
* There should be a function to check whether a BDD respects a given variable ordering
* There should be a function to `reduce` a BDD w.r.t. a specified variable ordering, yielding a ROBDD.

Desirable features:

* A constructor that compiles a ROBDD from a Boolean expression
* A DSL for constructing ROBDDs from Boolean functions
* A DSL 
