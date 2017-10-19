#pragma once
#include "../StableFluid3d.h"

using namespace System;
using namespace System::ComponentModel::Composition;
using namespace VVVV::PluginInterfaces::V2;
using namespace VVVV::Utils::Streams;

namespace VVVV
{
	namespace StableFluidNodes
	{
		/// \brief managed wrapper class for vvvv
		public value class Fluid3d
		{
		public:
			/// \brief unmanaged 3d fluid pointer
			StableFluid3d * Fluid;
		};

		[PluginInfo(Name = "StableFluid", Category = "StableFluid", Version = "3d", Author = "woei", Credits = "Jos Stam",
					Help = "Creates and holds a fluid simulation space")]
		public ref class CreateFluid3d : IPluginEvaluate
		{
		public:
			[Input("Resolution X", DefaultValue = 2)]
			IDiffSpread<int>^ FResX;
			[Input("Resolution Y", DefaultValue = 2)]
			IDiffSpread<int>^ FResY;
			[Input("Resolution Z", DefaultValue = 2)]
			IDiffSpread<int>^ FResZ;
			[Input("Width", DefaultValue = 1)]
			IDiffSpread<float>^ FWidth;
			[Input("Height", DefaultValue = 1)]
			IDiffSpread<float>^ FHeight;
			[Input("Depth", DefaultValue = 1)]
			IDiffSpread<float>^ FDepth;
			[Input("Reset")]
			ISpread<bool>^ FReset;

			[Output("Fluid")]
			ISpread<Fluid3d>^ FOutput;
			[Output("Changed")]
			ISpread<bool>^ FChanged;

			virtual void Evaluate(int spreadMax);

			~CreateFluid3d();
		};

		[PluginInfo(Name = "Evaluate", Category = "StableFluid", Version = "3d", Author = "woei", Credits = "Jos Stam",
					Help = "Drives the fluid simulation, giving access to global parameters")]
		public ref class EvaluateFluid3d : IPluginEvaluate
		{
		public:
			[Input("Fluid")]
			ISpread<Fluid3d>^ FInput;

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
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};
		
		[PluginInfo(Name = "SetVelocity", Category = "StableFluid", Version = "3d", Author = "woei",
					Help = "Sets the velocities of the entire fluid simulation space")]
		public ref class SetVelocityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("U")]
			IInStream<float>^ FU;
			[Input("V")]
			IInStream<float>^ FV;
			[Input("W")]
			IInStream<float>^ FW;
			[Input("Set")]
			ISpread<bool>^ FSet;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};
		
		[PluginInfo(Name = "AddVelocity", Category = "StableFluid", Version = "3d", Author = "woei",
			Help = "Adds velocities to the entire fluid simulation space")]
		public ref class AddVelocityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("U")]
			IInStream<float>^ FU;
			[Input("V")]
			IInStream<float>^ FV;
			[Input("W")]
			IInStream<float>^ FW;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "AddSphericalVelocity", Category = "StableFluid", Version = "3d", Author = "woei",
					Help = "Adds velocities around the given center to the fluid simulation space")]
		public ref class AddSphericalVelocityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("Center X")]
			ISpread<float>^ FX;
			[Input("Center Y")]
			ISpread<float>^ FY;
			[Input("Center Z")]
			ISpread<float>^ FZ;
			[Input("Scale X", DefaultValue = 1)]
			ISpread<float>^ FSX;
			[Input("Scale Y", DefaultValue = 1)]
			ISpread<float>^ FSY;
			[Input("Scale Z", DefaultValue = 1)]
			ISpread<float>^ FSZ;
			[Input("Radius", DefaultValue = 0.1)]
			ISpread<float>^ FRad;
			[Input("Ramp")]
			ISpread<float>^ FRamp;
			[Input("U")]
			ISpread<float>^ FU;
			[Input("V")]
			ISpread<float>^ FV;
			[Input("W")]
			ISpread<float>^ FW;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "AddLocalVelocity", Category = "StableFluid", Version = "3d", Author = "woei")]
		public ref class AddLocalVelocityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("U")]
			IInStream<float>^ FU;
			[Input("V")]
			IInStream<float>^ FV;
			[Input("W")]
			IInStream<float>^ FW;
			[Input("Position X")]
			IInStream<float>^ FX;
			[Input("Position Y")]
			IInStream<float>^ FY;
			[Input("Position Z")]
			IInStream<float>^ FZ;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "SetDensity", Category = "StableFluid", Version = "3d", Author = "woei",
					Help = "Sets the densities of the entire fluid simulation space")]
		public ref class SetDensityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("Density")]
			IInStream<float>^ FDens;
			[Input("Set")]
			ISpread<bool>^ FSet;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "AddDensity", Category = "StableFluid", Version = "3d", Author = "woei",
					Help = "Adds densities to the entire fluid simulation space")]
		public ref class AddDensityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("Density")]
			IInStream<float>^ FDens;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "AddSphericalDensity", Category = "StableFluid", Version = "3d", Author = "woei",
					Help = "Adds densities around the given center to the fluid simulation space")]
		public ref class AddSphericalDensityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("Center X")]
			ISpread<float>^ FX;
			[Input("Center Y")]
			ISpread<float>^ FY;
			[Input("Center Z")]
			ISpread<float>^ FZ;
			[Input("Scale X", DefaultValue = 1)]
			ISpread<float>^ FSX;
			[Input("Scale Y", DefaultValue = 1)]
			ISpread<float>^ FSY;
			[Input("Scale Z", DefaultValue = 1)]
			ISpread<float>^ FSZ;
			[Input("Radius", DefaultValue = 0.1)]
			ISpread<float>^ FRad;
			[Input("Ramp")]
			ISpread<float>^ FRamp;
			[Input("Density")]
			ISpread<float>^ FDens;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "AddLocalDensity", Category = "StableFluid", Version = "3d", Author = "woei")]
		public ref class AddLocalDensityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("Density")]
			IInStream<float>^ FDens;
			[Input("Position X")]
			IInStream<float>^ FX;
			[Input("Position Y")]
			IInStream<float>^ FY;
			[Input("Position Z")]
			IInStream<float>^ FZ;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "AddSpherical", Category = "StableFluid", Version = "3d", Author = "woei",
					Help = "Adds velocities and densities around the given center to the fluid simulation space")]
		public ref class AddSphericalFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("Center X")]
			ISpread<float>^ FX;
			[Input("Center Y")]
			ISpread<float>^ FY;
			[Input("Center Z")]
			ISpread<float>^ FZ;
			[Input("Scale X", DefaultValue = 1)]
			ISpread<float>^ FSX;
			[Input("Scale Y", DefaultValue = 1)]
			ISpread<float>^ FSY;
			[Input("Scale Z", DefaultValue = 1)]
			ISpread<float>^ FSZ;
			[Input("Radius", DefaultValue = 0.1)]
			ISpread<float>^ FRad;
			[Input("Ramp")]
			ISpread<float>^ FRamp;
			[Input("U")]
			ISpread<float>^ FU;
			[Input("V")]
			ISpread<float>^ FV;
			[Input("W")]
			ISpread<float>^ FW;
			[Input("Density")]
			ISpread<float>^ FDens;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "AddLocal", Category = "StableFluid", Version = "3d", Author = "woei")]
		public ref class AddLocalFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Input("U")]
			IInStream<float>^ FU;
			[Input("V")]
			IInStream<float>^ FV;
			[Input("W")]
			IInStream<float>^ FW;
			[Input("Density")]
			IInStream<float>^ FDens;
			[Input("Position X")]
			IInStream<float>^ FX;
			[Input("Position Y")]
			IInStream<float>^ FY;
			[Input("Position Z")]
			IInStream<float>^ FZ;
			[Input("Add")]
			ISpread<bool>^ FAdd;

			[Output("Output")]
			ISpread<Fluid3d>^ FOutput;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "GetVelocity", Category = "StableFluid", Version = "3d", Author = "woei",
					Help = "Returns the velocities of the enire fluid simulation space")]
		public ref class GetVelocityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Output("Fluid")]
			ISpread<Fluid3d>^ FOutput;
			[Output("U")]
			IOutStream<float>^ FU;
			[Output("V")]
			IOutStream<float>^ FV;
			[Output("W")]
			IOutStream<float>^ FW;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "GetDensity", Category = "StableFluid", Version = "3d", Author = "woei",
					Help = "Returns the densities of the enire fluid simulation space")]
		public ref class GetDensityFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Output("Fluid")]
			ISpread<Fluid3d>^ FOutput;
			[Output("Density")]
			IOutStream<float>^ FDens;

			virtual void Evaluate(int spreadMax);
		};

		[PluginInfo(Name = "Split", Category = "StableFluid", Version = "3d", Author = "woei",
					Help = "Returns the velocities and the densities of the enire fluid simulation space")]
		public ref class SplitFluid3d : IPluginEvaluate
		{
		public:
			[Input("Input")]
			ISpread<Fluid3d>^ FInput;

			[Output("U")]
			IOutStream<float>^ FU;
			[Output("V")]
			IOutStream<float>^ FV;
			[Output("W")]
			IOutStream<float>^ FW;
			[Output("Density")]
			IOutStream<float>^ FDens;

			virtual void Evaluate(int spreadMax);
		};
	}
}

