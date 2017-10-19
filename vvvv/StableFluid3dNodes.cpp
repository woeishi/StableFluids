#include "../Stdafx.h"
#include "StableFluid3dNodes.h"
#include "../StableFluid3d.h"

using namespace System::ComponentModel::Composition;
using namespace VVVV::PluginInterfaces::V2;

namespace VVVV
{
	namespace StableFluidNodes
	{
		void CreateFluid3d::Evaluate(int spreadMax)
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
				if (FResX->IsChanged || FResY->IsChanged || FResZ->IsChanged
					|| FWidth->IsChanged || FHeight->IsChanged || FDepth->IsChanged || FReset[0])
				{
					FOutput->SliceCount = 1;
					if (FOutput[0].Fluid)
					{
						delete FOutput[0].Fluid;
					}
					ISpread<int>^ resX = FResX;
					ISpread<int>^ resY = FResY;
					ISpread<int>^ resZ = FResZ;
					ISpread<float>^ width = FWidth;
					ISpread<float>^ height = FHeight;
					ISpread<float>^ depth = FDepth;

					Fluid3d f;
					f.Fluid = new StableFluid3d(resX[0], resY[0], resZ[0], width[0], height[0], depth[0]);
					FOutput[0] = f;
					FChanged[0] = true;
				}
				else
					FChanged[0] = false;
			}
		}
		CreateFluid3d::~CreateFluid3d()
		{
			if (FOutput->SliceCount > 0 && FOutput[0].Fluid)
			{
				delete FOutput[0].Fluid;
			}
		}
		
		void EvaluateFluid3d::Evaluate(int spreadMax)
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

		void SetVelocityFluid3d::Evaluate(int spreadMax)
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
					auto _w = StreamUtils::GetCyclicReader(FW);

					float *u0 = f->U();
					float *v0 = f->V();
					float *w0 = f->W();
					for (int l = 0; l < f->Count(); u0++, v0++, w0++, l++)
					{
						*u0 += _u->Read(1);
						*v0 += _v->Read(1);
						*w0 += _w->Read(1);
					}
					delete _u;
					delete _v;
					delete _w;
				}
			}
		}
	
		void AddVelocityFluid3d::Evaluate(int spreadMax)
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
					auto _w = StreamUtils::GetCyclicReader(FW);

					float *u0 = f->NewU();
					float *v0 = f->NewV();
					float *w0 = f->NewW();
					for (int l = 0; l < f->Count(); u0++, v0++, w0++, l++)
					{
						if (!_u->Eos)
							*u0 += _u->Read(1);
						if (!_v->Eos)
							*v0 += _v->Read(1);
						if (!_w->Eos)
							*w0 += _w->Read(1);
					}
					delete _u;
					delete _v;
					delete _w;
				}
			}
		}
		
		void AddSphericalVelocityFluid3d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					float * u = FOutput[0].Fluid->NewU();
					float * v = FOutput[0].Fluid->NewV();
					float * w = FOutput[0].Fluid->NewW();
					float * target[] = { u, v, w };
					for (int i = 0; i < spreadMax; i++)
					{
						float val[] = { FU[i], FV[i], FW[i] };
						FOutput[0].Fluid->AddSpherical(FX[i], FY[i], FZ[i], FSX[i], FSY[i], FSZ[i], FRad[i], FRamp[i], 3, val, target);
					}
				}
			}
		}
		
		void AddLocalVelocityFluid3d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					auto f = FOutput[0].Fluid;
					auto _x = FX->GetReader(), _y = FY->GetReader(), _z = FZ->GetReader();
					auto _u = FU->GetReader(), _v = FV->GetReader(), _w = FW->GetReader();
					while ((!_x->Eos) && (!_y->Eos) && (!_z->Eos))
					{
						float x = _x->Read(1), y = _y->Read(1), z = _z->Read(1);
						if (!_u->Eos)
							f->TrilinearAdd(x, y, z, _u->Read(1), f->NewU());
						if (!_v->Eos)
							f->TrilinearAdd(x, y, z, _v->Read(1), f->NewV());
						if (!_w->Eos)
							f->TrilinearAdd(x, y, z, _w->Read(1), f->NewW());
					}
					delete _x, _y, _z, _u, _v, _w;
				}
			}
		}
		
		void SetDensityFluid3d::Evaluate(int spreadMax)
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
	
		void AddDensityFluid3d::Evaluate(int spreadMax)
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

		void AddSphericalDensityFluid3d::Evaluate(int spreadMax)
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
						float val[] = { FDens[i] };
						FOutput[0].Fluid->AddSpherical(FX[i], FY[i], FZ[i], FSX[i], FSY[i], FSZ[i], FRad[i], FRamp[i], 1, val, target);
					}
				}
			}
		}

		void AddLocalDensityFluid3d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					auto f = FOutput[0].Fluid;
					auto _x = FX->GetReader(), _y = FY->GetReader(), _z = FZ->GetReader();
					auto _d = FDens->GetReader();
					while ((!_x->Eos) && (!_y->Eos) && (!_z->Eos) && (!_d->Eos))
					{
						f->TrilinearAdd(_x->Read(1), _y->Read(1), _z->Read(1), _d->Read(1), f->NewDensity());
					}
					delete _x, _y, _z, _d;
				}
			}
		}

		void AddSphericalFluid3d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					float * u = FOutput[0].Fluid->NewU();
					float * v = FOutput[0].Fluid->NewV();
					float * w = FOutput[0].Fluid->NewW();
					float * d = FOutput[0].Fluid->NewDensity();
					float * target[] = { u, v, w, d };
					for (int i = 0; i < spreadMax; i++)
					{
						float val[] = { FU[i], FV[i], FW[i], FDens[i] };
						FOutput[0].Fluid->AddSpherical(FX[i], FY[i], FZ[i], FSX[i], FSY[i], FSZ[i], FRad[i], FRamp[i], 4, val, target);
					}
				}
			}
		}

		void AddLocalFluid3d::Evaluate(int spreadMax)
		{
			FOutput->SliceCount = FInput->SliceCount;
			if (FInput->SliceCount > 0)
			{
				FOutput[0] = FInput[0];
				if (FOutput[0].Fluid && spreadMax > 0 && FAdd[0])
				{
					auto f = FOutput[0].Fluid;
					auto _x = FX->GetReader(), _y = FY->GetReader(), _z = FZ->GetReader();
					auto _u = FU->GetReader(), _v = FV->GetReader(), _w = FW->GetReader(), _d = FDens->GetReader();
					while ((!_x->Eos) && (!_y->Eos) && (!_z->Eos))
					{
						float x = _x->Read(1), y = _y->Read(1), z = _z->Read(1);
						if (!_u->Eos)
							f->TrilinearAdd(x, y, z, _u->Read(1), f->NewU());
						if (!_v->Eos)
							f->TrilinearAdd(x, y, z, _v->Read(1), f->NewV());
						if (!_w->Eos)
							f->TrilinearAdd(x, y, z, _w->Read(1), f->NewW());
						if (!_d->Eos)
							f->TrilinearAdd(x, y, z, _d->Read(1), f->NewDensity());
					}
					delete _x, _y, _z, _u, _v, _w, _d;
				}
			}
		}

		void GetVelocityFluid3d::Evaluate(int spreadMax)
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
					((IOutStream^)FW)->Length = f->Count();
					auto _u = FU->GetWriter(), _v = FV->GetWriter(), _w = FW->GetWriter();
					float *outU = f->U(), *outV = f->V(), *outW = f->W();
					for (int i = 0; i < f->Count(); outU++, outV++, outW++, i++)
					{
						_u->Write(*outU, 1);
						_v->Write(*outV, 1);
						_w->Write(*outW, 1);
					}
					delete _u, _v, _w;
				}
				else
				{
					((IOutStream^)FU)->Length = 0;
					((IOutStream^)FV)->Length = 0;
					((IOutStream^)FW)->Length = 0;
				}
			}
		}

		void GetDensityFluid3d::Evaluate(int spreadMax)
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

		void SplitFluid3d::Evaluate(int spreadMax)
		{
			if (FInput->SliceCount >0 && FInput[0].Fluid)
			{
				auto f = FInput[0].Fluid;
				((IOutStream^)FU)->Length = f->Count();
				((IOutStream^)FV)->Length = f->Count();
				((IOutStream^)FW)->Length = f->Count();
				((IOutStream^)FDens)->Length = f->Count();
				auto _u = FU->GetWriter(), _v = FV->GetWriter(), _w = FW->GetWriter();
				auto _d = FDens->GetWriter();

				float *outU = f->U(), *outV = f->V(), *outW = f->W();
				float *outD = f->Density();
				for (int i = 0; i<f->Count(); outU++, outV++, outW++, outD++, i++)
				{
					_u->Write(*outU, 1);
					_v->Write(*outV, 1);
					_w->Write(*outW, 1);
					_d->Write(*outD, 1);
				}
				delete _u, _v, _w;
				delete _d;
			}
			else
			{
				((IOutStream^)FU)->Length = 0;
				((IOutStream^)FV)->Length = 0;
				((IOutStream^)FW)->Length = 0;
				((IOutStream^)FDens)->Length = 0;
			}
		}	
	}
}