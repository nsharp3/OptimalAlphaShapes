OptimalAlphaShapes
==================

Python implementation of optimal trajectory calculation in unsteady flows using alpha-shapes. Currently a work in progress.

**Authors:** Nicholas Sharp (nsharp3@vt.edu) and Shane Ross (sdross@vt.edu)

###Overview
Uses alpha-shapes to mesh a 3D front-propagation search. See paper for more details (TODO: paper URL). Solves two main problems:

* Finding time-optimal paths in 3D current fields
* Finding optimal paths in 2d current fields when the optimality metric includes both time and energy

Note in both cases, this method is valid on time-varying fields, as well as the weakly-propelled case where the currents are significantly stronger than the vehicle's propulsion.

To the authors knowledge, this is the first published method for solving either of these problems.

###Current questions to be resolved
* Is compounding numerical error introduced when the controls are integrated from the ODEs instead of re-calculated geometrically from the mesh?
* Classification of facets in the alpha shape from CGAL. Sometime extra interior facets are retained, is there a way to use the classifications to avoid this?
* Completion testing. Partly because of the facet classification issue, completition is tested simply by nearness to the hull. There are many other ways to do this.
* Overall, this code is not tested as a system like I wish it was, basically because there is nothing to compare against. I welcome suggestion for known or verifiable results to use for testing.

###Running
Note that .cpp code which interfaces with CGAL needs to be compiled and linked against CGAL. This can be messy, you're your own for that. Look in to CGAL CMake scripts and how they work.

This is a RESEARCH code. It is meant as a proptype implementation of a new idea. It likely contains bugs. Do not trust this code for anything unless you have verfied its functionality yourself.
