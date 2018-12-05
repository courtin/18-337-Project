# 18-337-Project

#MIGP.jl

A Julia package that uses GPkit to solve Mixed-Integer Geometric Programs. 

Usage: Problems are formulated as GPkit-compatible optimization problems (see evtol.py and mico_test.py for examples). Integer varables are declared by adding "INT_" to the variable name in the docstring declaration. For example, declaring "INT_myvar" tells MIGP.jl to treat myvar as an integer variable.  

Currently, integer variable declarations must be in the top level model.  

**Julia Packages:** LightGraphs, TikzGraphs, TikzPictures, PyCall, JuMP, Mosek, Distributed

**Python Packages:** gpkit, gpfit, gplibrary, mosek/cvxopt

Installation instructions for gpkit and Mosek can be found here: https://gpkit.readthedocs.io/en/latest/

This package is solver agnostic; Mosek requires a (free academic) license.  cvxopt is open source and will also work (albeit more slowly). 

gpfit and gplibrary are also required, and can be found at https://github.com/convexengineering.

