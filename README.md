# ROBDD.jl

A Julia implementation of (Reduced Ordered) [Binary Decision Diagrams](https://en.wikipedia.org/wiki/Binary_decision_diagram).

> "...one of the only really fundamental data structures that came out in the last twenty-five years"
>
> -- Donald Knuth. [_Fun with Binary Decision Diagrams_](https://www.youtube.com/watch?v=SQE21efsf7Y). 2008


## Overview

My design choices are informed by this document:

[An Introduction to Binary Decision Diagrams](https://www.cmi.ac.in/~madhavan/courses/verification-2011/andersen-bdd.pdf). Henrik Reif Andersen. 1999.

(However, my basic ROBDD construction algorithm builds up the ROBDD from boolean operations, via DFS.
The `build` algorithm mentioned in that document includes substitution operations which are costly to implement.)

Some basic requirements:

* Build an ROBDD from a Boolean expression
* Perform boolean operations on ROBDDs
* Restrict ROBDD
* Sample an arbitrary model from an ROBDD
* Enumerate models of an ROBDD

Desirable features:

* A DSL for constructing Boolean expressions
* Count models of an ROBDD
* Variable reordering
* Simplify an ROBDD w.r.t. a "domain of interest"
* Existential quantification 

## Installation

## Usage
