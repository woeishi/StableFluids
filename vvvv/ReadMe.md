#StableFluidsNodes

C++/CLI VVVV wrapper classes/nodes

for 2d and 3d version alike the value class FluidNd just holds a pointer to the unmanaged StableFluid instance.
one managed class instance only (CreateFluidNd) can create and delete it.

uses the NugetDummy project to pull the .NET dependencies

Nodes
---
* __StableFluid__: Creates and holds a fluid simulation space
* __Evaluate__: Drives the fluid simulation, giving access to global parameters
* __SetVelocity / SetDensity__: Sets the velocities/densities of the entire fluid simulation space
* __AddVelocity / AddDensity__: Adds velocities/densities to the entire fluid simulation space
* __AddRadialVelocity / AddRadialDensity__: Adds velocities/densities around the given center to the fluid simulation space (think of it as an emitter)
* __AddRadial__: Adds velocities and densities around the given center to the fluid simulation space
* __GetVelocity / GetDensity__: Returns the velocities/densities of the enire fluid simulation space
* __Split__: Returns the velocities and densities of the enire fluid simulation space

