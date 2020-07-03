# BDD.jl

A Julia implementation of [Binary Decision Diagrams](https://en.wikipedia.org/wiki/Binary_decision_diagram).

> "...one of the only really fundamental data structures that came out in the last twenty-five years"
>
> -- Donald Knuth. [_Fun with Binary Decision Diagrams_](https://www.youtube.com/watch?v=SQE21efsf7Y). 2008


Some basic requirements:

* User should be able to manually construct a BDD by specifying decision nodes
* User should be able to `reduce` a BDD w.r.t. a specified variable ordering

Some desirable features:

* A constructor that compiles a ROBDD from a Boolean expression
* A DSL that 

