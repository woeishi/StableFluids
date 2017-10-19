#include "../Stdafx.h"
#include "StableFluid2dNodes.h"
#include "../StableFluid2d.h"

using namespace System::ComponentModel::Composition;
using namespace VVVV::PluginInterfaces::V2;


namespace VVVV
{
	namespace StableFluidNodes
	{
		void CreateFluid2d::Evaluate(int spreadMax)
		{
			if (spreadMax == 0)
			{
				if (FOutput[0].Fluid)
				{
					delete FOutput[0].Fluid;
				}
				FOutput->SliceCount = 0;
				FChanged->SliceCount = 0;
			}
			else
			{
				FChanged->SliceCount = 1;
				if (FResX->IsChanged || FResY->IsChanged
					|| FWidth->IsChanged || FHeight->IsChanged || FReset[0])
				{
					FOutput->SliceCount = 1;
					if (FOutput[0].Fluid)
					{
						delete FOutput[0].Fluid;
					}
					ISpread<int>^ resX = FResX;
					ISpread<int>^ resY = FResY;
					ISpread<float>^ width = FWidth;
					ISpread<float>^ height = FHeight;

					Fluid2d f;
					f.Fluid = new StableFluid2d((int)resX[0], (int)resY[0], (float)width[0], (float)height[0]);
					FOutput[0] = f;
					FChanged[0] = true;
				}
				else
					FChanged[0] = false;
			}
		}
		CreateFluid2d::~CreateFluid2d()
		{
			if (FOutput->SliceCount > 0 && FOutput[0].Fluid)
			{
				delete FOutput[0].Fluid;
			}
		}
		
		void EvaluateFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FRun[0])
				{
					FOutput[0].Fluid->SolverIterations = FSolvIter[0];
					if (FMulVel[0] != 1)
						FOutput[0].Fluid->MultiplyVelocity(FMulVel[0]);
					if (FMulDens[0] != 1)
						FOutput[0].Fluid->MultiplyDensity(FMulDens[0]);
					FOutput[0].Fluid->Update(FVisc[0], FDiff[0], FDt[0]);
				}
			}
		}

		void SetVelocityFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FSet[0])
				{
					auto f = FOutput[0].Fluid;
					auto _u = StreamUtils::GetCyclicReader(FU);
					auto _v = StreamUtils::GetCyclicReader(FV);

					float *u0 = f->U();
					float *v0 = f->V();
					for (int l = 0; l < f->Count(); u0++, v0++, l++)
					{
						*u0 += _u->Read(1);
						*v0 += _v->Read(1);
					}
					delete _u, _v;
				}
			}
		}

		void AddVelocityFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					auto f = FOutput[0].Fluid;
					auto _u = StreamUtils::GetCyclicReader(FU);
					auto _v = StreamUtils::GetCyclicReader(FV);

					float* u0 = f->NewU();
					float* v0 = f->NewV();
					for (int l = 0; l < f->Count(); u0++, v0++, l++)
					{
						*u0 += _u->Read(1);
						*v0 += _v->Read(1);
					}
					delete _u, _v;
				}
			}
		}

		void AddRadialVelocityFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					float * u = FOutput[0].Fluid->NewU();
					float * v = FOutput[0].Fluid->NewV();
					float * target[] = { u, v };
					for (int i = 0; i < spreadMax; i++)
					{
						float val[] = { FU[i], FV[i] };
						
						FOutput[0].Fluid->AddRadial(FX[i], FY[i], FSX[i], FSY[i], FRad[i], FRamp[i], 2, val, target);
					}
				}
			}
		}

		void SetDensityFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FSet[0])
				{
					auto f = FOutput[0].Fluid;
					auto _d = StreamUtils::GetCyclicReader(FDens);
					int i = 0;
					float* dens0 = f->Density();
					while (i < f->Count())
					{
						*dens0 += _d->Read(1);
						dens0++;
						i++;
					}
					delete _d;
				}
			}
		}

		void AddDensityFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					auto f = FOutput[0].Fluid;
					auto _d = StreamUtils::GetCyclicReader(FDens);
					int i = 0;
					float* dens0 = f->NewDensity();
					while (i < f->Count())
					{
						*dens0 += _d->Read(1);
						dens0++;
						i++;
					}
					delete _d;
				}
			}
		}

		void AddRadialDensityFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					float * d = FOutput[0].Fluid->NewDensity();
					float * target[] = { d };
					for (int i = 0; i < spreadMax; i++)
					{
						float val[]= { FDens[i] };
						FOutput[0].Fluid->AddRadial(FX[i], FY[i], FSX[i], FSY[i], FRad[i], FRamp[i], 1, val, target);
					}
				}
			}
		}

		void AddRadialFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					float * u = FOutput[0].Fluid->NewU();
					float * v = FOutput[0].Fluid->NewV();
					float * d = FOutput[0].Fluid->NewDensity();
					float * target[] = { u, v, d };
					for (int i = 0; i < spreadMax; i++)
					{
						float val[] = { FU[i], FV[i], FDens[i] };
						FOutput[0].Fluid->AddRadial(FX[i], FY[i], FSX[i], FSY[i], FRad[i], FRamp[i], 3, val, target);
					}
				}
			}
		}

		void GetVelocityFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid)
				{
					auto f = FOutput[0].Fluid;
					((IOutStream^)FU)->Length = f->Count();
					((IOutStream^)FV)->Length = f->Count();
					auto _u = FU->GetWriter(), _v = FV->GetWriter();
					float *outU = f->U(), *outV = f->V();
					for (int i = 0; i < f->Count(); outU++, outV++, i++)
					{
						_u->Write(*outU, 1);
						_v->Write(*outV, 1);
					}
					delete _u, _v;
				}
				else
				{
					((IOutStream^)FU)->Length = 0;
					((IOutStream^)FV)->Length = 0;
				}
			}
		}

		void GetDensityFluid2d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid)
				{
					auto f = FOutput[0].Fluid;
					((IOutStream^)FDens)->Length = f->Count();
					auto _d = FDens->GetWriter();
					float* outD = f->Density();
					for (int i = 0; i < f->Count(); outD++, i++)
						_d->Write(*outD, 1);
					delete _d;
				}
				else
				{
					((IOutStream^)FDens)->Length = 0;
				}
			}
		}
		
		void SplitFluid2d::Evaluate(int spreadMax)
		{
			if (FInput->SliceCount > 0 && FInput[0].Fluid)
			{
				auto f = FInput[0].Fluid;
				((IOutStream^)FU)->Length = f->Count();
				((IOutStream^)FV)->Length = f->Count();
				((IOutStream^)FDens)->Length = f->Count();
				auto u = FU->GetWriter();
				auto v = FV->GetWriter();
				auto d = FDens->GetWriter();

				float* outU = f->U();
				float* outV = f->V();
				float* outD = f->Density();
				for (int i = 0; i < f->Count(); outU++, outV++, outD++, i++)
				{
					u->Write(*outU, 1);
					v->Write(*outV, 1);
					d->Write(*outD, 1);
				}
				delete u, v;
				delete d;
			}
		}
	}
}