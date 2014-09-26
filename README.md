OptimalAlphaShapes
==================

Python implementation of optimal trajectory calculation in unsteady flows using alpha-shapes.

**Authors:** Nicholas Sharp (nsharp3@vt.edu) and Shane Ross (sdross@vt.edu)

###Overview
Uses alpha-shapes to mesh a 3D front-propagation search. See paper for more details (submitted to ACC 2015). Solves two main problems:

* Finding time-optimal paths in 3D current fields
* Finding optimal paths in 2d current fields when the optimality metric includes both time and energy

(Only the first is addressed in the ACC submission, the latter is the focus of upcoming work)

Note in both cases, this method is valid on time-varying fields, as well as the weakly-propelled case where the currents are significantly stronger than the vehicle's propulsion.

To the authors' knowledge, this is the first published method for solving either of these problems

###Running
Note that the CGAL_Alpha_Wrapper.cpp code which interfaces with CGAL needs to be compiled and linked against CGAL. This can be difficult, and the authors do not provide support. Look in to CGAL CMake scripts and how they work. A binary is provided which was compiled on 64-bit linux, there is a small chance this may work for you but it is very unlikely.


*This is a RESEARCH code. It is meant as a prototype implementation of a new idea. It may contains bugs. Do not trust this code for anything important unless you have verfied its functionality yourself.*
