StableFluids
===

basic **real-time fluids** with SIMD optimization
direct implementation of [Jos Stam's Stable Fluid algorithm](http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf)

> This project quite strictly reproduces the contents of the paper since it mainly served as a __reference__ for evaluating other algorithms and implementation attempts.
> The only deviations are the additional vorticity confinement routines, eliminating the restriction to a square and evenly spaced base grid, and using Jacobi instead of Gauss-Seidel solver in order to benefit from SIMD instructions.

features
---

*  2d and 3d versions
*  independent simulation space resolution and grid spacing per dimension
*  vorticity confinement
*  viscosity, diffusion and timestep parameter
*  variable solver iteration count with two extra parameters to manually balance the simulation even if the count is set too low
*  VVVV wrapper nodes 

structure
---
* unmanaged c/c++ fluid classes (must be unmanaged for SSE/AVX)
* CLI wrapper classes for VVVV nodes (must be managed for VVVV)
* a C# project to pull the NuGet packages required for VVVV

so while this is a C++/CLI Visual Studio project, StableFluid2d.cpp and StableFluid3d.cpp are set to build without Common Runtime Language support in release mode to enable AVX instructions.

###sidenote
the structure and the code might look a bit funny because of reasons:

* i usually write C#
* first time seriously doing a c++ project from scratch
* first time peaking into SIMD instructions
* looking for a nice way combining VVVVs NET. plugininterface with optimized c++

while the datastructure looks like it could just benefit from automatic modern compiler optimizations when having the the floats nicely aligned, it's actually the border conditions, which break up the alignment and therefore require manual vectorization in a lot of the base routines.

