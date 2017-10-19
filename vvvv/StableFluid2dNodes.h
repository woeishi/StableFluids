#pragma once
#include "../StableFluid2d.h"

using namespace System;
using namespace System::ComponentModel::Composition;
using namespace VVVV::PluginInterfaces::V2;
using namespace VVVV::Utils::Streams;

namespace VVVV
{
	namespace StableFluidNodes
	{
		/// \brief managed wrapper class for vvvv
		public value class Fluid2d
		{
		public:
			/// \brief unmanaged 2d fluid pointer
			StableFluid2d * Fluid;
		};

		[PluginInfo(Name = "StableFluid", Category = "StableFluid", Version = "2d", Author = "woei", Credits = "Jos Stam",
					Help = "Creates and holds a fluid simulation space")]
		public ref class CreateFluid2d : IPluginEvaluate
		{
		public:
			[Input("Resolution X", DefaultValue = 2)]
			IDiffSpread<int>^ FResX;
			[Input("Resolution Y", DefaultValue = 2)]
			IDiffSpread<int>^ FResY;
			[Input("Width", DefaultValue = 1)]
			IDiffSpread<float>^ FWidth;
			[Input("Height", DefaultValue = 1)]
			IDiffSpread<float>^ FHeight;
			[Input("Reset")]
			ISpread<bool>^ FReset;

			[Output("Fluid")]
			ISpread<Fluid2d>^ FOutput;
			[Output("Changed")]
			ISpread<bool>^ FChanged;

			virtual void Evaluate(int spreadMax);

			~CreateFluid2d();
		};

		[PluginInfo(Name = "Evaluate", Category = "StableFluid", Version = "2d", Author = "woei", Credits = "Jos Stam",
					Help = "Drives the fluid simulation, giving access to global parameters")]
		public ref class EvaluateFluid2d : IPluginEvaluate
		{
		public:
			[Input("Fluid")]
			ISpread<Fluid2d>^ FInput;

			[Input("Viscosity")]
			ISpread<float>^ FVisc;
			[Input("Multiply Velocity", DefaultValue = 1, Visibility = PinVisibility::OnlyInspector)]
			ISpread<float>^ FMulVel;
			[Input("Diffusion")]
			ISpread<float>^ FDiff;
			[Input("Multiply Density", DefaultValue = 1, Visibility = PinVisibility::OnlyInspector)]
			ISpread<float>^ FMulDens;
			[Input("Delta Time", DefaultValue = 0.016666)]
			ISpread<float>^ FDt;
			[Input("Solver Iterations", DefaultValue = 20, Visibility = PinVisibility::OnlyInspector)]
			ISpread<int>^ FSolvIter;
			[Input("Update")]
			ISpread<bool>^ FRun;

			[Output("Fluid")]
			ISpread<Fluid2d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "SetVelocity", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Sets the velocities of the entire fluid simulation space")]
		public ref class SetVelocityFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Input("U")]
			IInStream<float>^ FU;
			[Input("V")]
			IInStream<float>^ FV;
			[Input("Set")]
			ISpread<bool>^ FSet;

			[Output("Output")]
			ISpread<Fluid2d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};
		
		[PluginInfo(Name = "AddVelocity", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Adds velocities to the entire fluid simulation space")]
		public ref class AddVelocityFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Input("U")]
			IInStream<float>^ FU;
			[Input("V")]
			IInStream<float>^ FV;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid2d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "AddRadialVelocity", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Adds velocities around the given center to the fluid simulation space")]
		public ref class AddRadialVelocityFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Input("Center X")]
			ISpread<float>^ FX;
			[Input("Center Y")]
			ISpread<float>^ FY;
			[Input("Scale X", DefaultValue = 1)]
			ISpread<float>^ FSX;
			[Input("Scale Y", DefaultValue = 1)]
			ISpread<float>^ FSY;
			[Input("Radius", DefaultValue = 0.1)]
			ISpread<float>^ FRad;
			[Input("Ramp")]
			ISpread<float>^ FRamp;
			[Input("U")]
			ISpread<float>^ FU;
			[Input("V")]
			ISpread<float>^ FV;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid2d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};
		
		[PluginInfo(Name = "SetDensity", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Sets the densities of the entire fluid simulation space")]
		public ref class SetDensityFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Input("Density")]
			IInStream<float>^ FDens;
			[Input("Set")]
			ISpread<bool>^ FSet;

			[Output("Output")]
			ISpread<Fluid2d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};
		
		[PluginInfo(Name = "AddDensity", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Adds densities to the entire fluid simulation space")]
		public ref class AddDensityFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Input("Density")]
			IInStream<float>^ FDens;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid2d>^ FOutput;

			virtual void Evaluate(int spreadMax);
			
		};
		
		[PluginInfo(Name = "AddRadialDensity", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Adds densities around the given center to the fluid simulation space")]
		public ref class AddRadialDensityFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Input("Center X")]
			ISpread<float>^ FX;
			[Input("Center Y")]
			ISpread<float>^ FY;
			[Input("Scale X", DefaultValue = 1)]
			ISpread<float>^ FSX;
			[Input("Scale Y", DefaultValue = 1)]
			ISpread<float>^ FSY;
			[Input("Radius", DefaultValue = 0.1)]
			ISpread<float>^ FRad;
			[Input("Ramp")]
			ISpread<float>^ FRamp;
			[Input("Density")]
			ISpread<float>^ FDens;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid2d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "AddRadial", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Adds velocities and densities around the given center to the fluid simulation space")]
		public ref class AddRadialFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Input("Center X")]
			ISpread<float>^ FX;
			[Input("Center Y")]
			ISpread<float>^ FY;
			[Input("Scale X", DefaultValue = 1)]
			ISpread<float>^ FSX;
			[Input("Scale Y", DefaultValue = 1)]
			ISpread<float>^ FSY;
			[Input("Radius", DefaultValue = 0.1)]
			ISpread<float>^ FRad;
			[Input("Ramp")]
			ISpread<float>^ FRamp;
			[Input("U")]
			ISpread<float>^ FU;
			[Input("V")]
			ISpread<float>^ FV;
			[Input("Density")]
			ISpread<float>^ FDens;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid2d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "GetVelocity", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Returns the velocities of the enire fluid simulation space")]
		public ref class GetVelocityFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Output("Fluid")]
			ISpread<Fluid2d>^ FOutput;
			[Output("U")]
			IOutStream<float>^ FU;
			[Output("V")]
			IOutStream<float>^ FV;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "GetDensity", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Returns the densities of the enire fluid simulation space")]
		public ref class GetDensityFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Output("Fluid")]
			ISpread<Fluid2d>^ FOutput;
			[Output("Density")]
			IOutStream<float>^ FDens;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "Split", Category = "StableFluid", Version = "2d", Author = "woei",
					Help = "Returns the velocities and the densities of the enire fluid simulation space")]
		public ref class SplitFluid2d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid2d>^ FInput;

			[Output("U")]
			IOutStream<float>^ FU;
			[Output("V")]
			IOutStream<float>^ FV;
			[Output("Density")]
			IOutStream<float>^ FDens;

			virtual void Evaluate(int spreadMax);
		};
	}
}
