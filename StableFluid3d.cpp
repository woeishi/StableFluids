#include "stdafx.h"
#include "StableFluid3d.h"

#include <algorithm>

StableFluid3d::StableFluid3d(int countX, int countY, int countZ, float width, float height, float depth)
{
	length = countX * countY * countZ;
	resX = countX;
	nX = countX - 2;
	resY = countY;
	nY = countY - 2;
	resZ = countZ;
	nZ = countZ - 2;

	rangeX = width;
	hx = width / countX;
	hxp = pow(hx, -2.0f);
	maxX = (width - hx) / 2.0f;

	rangeY = height;
	hy = height / countY;
	hyp = pow(hy, -2.0f);
	maxY = (height - hy) / 2.0f;

	rangeZ = depth;
	hz = depth / countZ;
	hzp = pow(hz, -2.0f);
	maxZ = (depth - hz) / 2.0f;

	_u = new float[length]();
	_v = new float[length]();
	_w = new float[length]();
	_u0 = new float[length]();
	_v0 = new float[length]();
	_w0 = new float[length]();
	_cu = new float[length]();
	_cv = new float[length]();
	_cw = new float[length]();

	_dens = new float[length]();
	_dens0 = new float[length]();

	_scratch = new float[length]();

	// SSE/AVX
	yf = resX;
	zf = resX*resY;
	auto dr = std::div(nX, 4); //1-nX inclusive
	sseX = dr.quot * 4 + 1; //1+ because inclusive

	dr = std::div(nX, 8);
	avxX = dr.quot * 8 + 1;
}

StableFluid3d::~StableFluid3d()
{
	delete[] _v, _u0, _cu;
	delete[] _u, _v0, _cv;
	delete[] _w, _w0, _cw;

	delete[] _dens, _dens0;
	delete[] _scratch;
}

void StableFluid3d::ResetInputs()
{
	memset(_dens0, 0, length * sizeof(float));
	memset(_u0, 0, length * sizeof(float));
	memset(_v0, 0, length * sizeof(float));
	memset(_w0, 0, length * sizeof(float));
}

void StableFluid3d::Update(const float visc, const float diff, const float dt)
{
	//velocity
	AddSource(_u, _u0, dt);
	AddSource(_v, _v0, dt);
	AddSource(_w, _w0, dt);

	CurlU();
	CurlV();
	CurlW();

	VorticityConfinement(_u0, _v0, _w0);
	AddSource(_u, _u0, dt);
	AddSource(_v, _v0, dt);
	AddSource(_w, _w0, dt);

	Diffuse(1, _u0, _u, visc, dt);
	Diffuse(2, _v0, _v, visc, dt);
	Diffuse(3, _w0, _w, visc, dt);

	Project(_u0, _v0, _w0, _u, _v);

	Advect(1, _u, _u0, _u0, _v0, _w0, dt);
	Advect(2, _v, _v0, _u0, _v0, _w0, dt);
	Advect(3, _w, _w0, _u0, _v0, _w0, dt);
	Project(_u, _v, _w, _u0, _v0);

	//density
	AddSource(_dens, _dens0, dt);
	Diffuse(0, _dens0, _dens, diff, dt);
	Advect(0, _dens, _dens0, _u, _v, _w, dt);

	ResetInputs();
}

void StableFluid3d::MultiplyVelocity(const float rate)
{
	for (int i = 0; i < length; i++)
	{
		_u[i] *= rate;
		_v[i] *= rate;
		_w[i] *= rate;
	}
}

void StableFluid3d::MultiplyDensity(const float rate)
{
	for (int i = 0; i < length; i++)
		_dens[i] *= rate;
}

void StableFluid3d::AddSpherical(float cx, float cy, float cz, float sx, float sy, float sz, float r, float io, int paramSize, const float * value, float ** target)
{
	auto dr = div(resX, 8);
	int max8 = dr.quot * 8;
	float off[8]{ 0,1,2,3,4,5,6,7 };
	int i = 0;
	float px, py, pz, l, mask, d, d0;

	float zPos = -rangeZ * 0.5f + hz * 0.5f;
	for (int z = 0; z < resZ; z++)
	{
		pz = (cz - zPos) / sz;
		float yPos = -rangeY * 0.5f + hy * 0.5f;
		for (int y = 0; y < resY; y++)
		{
			py = (cy - yPos) / sy;
			float xPos = -rangeX * 0.5f + hx * 0.5f;
			for (int x = 0; x < max8; x += 8, i += 8)
			{
				__m256 _xPos = _mm256_add_ps(_mm256_set1_ps(xPos), _mm256_mul_ps(_mm256_loadu_ps(off), _mm256_set1_ps(hx)));
				__m256 _px = _mm256_div_ps(_mm256_sub_ps(_mm256_set1_ps(cx), _xPos), _mm256_set1_ps(sx));
				_px = _mm256_mul_ps(_px, _px);
				__m256 _py = _mm256_mul_ps(_mm256_set1_ps(py), _mm256_set1_ps(py));
				__m256 _pz = _mm256_mul_ps(_mm256_set1_ps(pz), _mm256_set1_ps(pz));
				__m256 _l = _mm256_sqrt_ps(_mm256_add_ps(_px, _mm256_add_ps(_py, _pz)));
				_l = _mm256_div_ps(_mm256_max_ps(_mm256_sub_ps(_mm256_set1_ps(r), _l), _mm256_set1_ps(0)), _mm256_set1_ps(r));
				__m256 _mask = _mm256_cmp_ps(_l, _mm256_set1_ps(0), _CMP_GT_OQ);

				for (int j = 0; j < paramSize; j++)
				{
					__m256 r1 = _mm256_mul_ps(_l, _mm256_mul_ps(_mm256_set1_ps(value[j]), _mm256_set1_ps(1 - io)));
					__m256 r2 = _mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps(1), _l), _mm256_mul_ps(_mm256_set1_ps(value[j]), _mm256_set1_ps(io)));
					r2 = _mm256_add_ps(r1, _mm256_and_ps(r2, _mask));
					_mm256_storeu_ps(target[j] + i, _mm256_add_ps(_mm256_loadu_ps(target[j] + i), r2));
				}
				xPos += hx * 8;
			}
			for (int x = max8; x < resX; x++, i++)
			{
				px = (cx - xPos) / sx;
				l = sqrt(px*px + py*py + pz*pz);
				l = fmax(r - l, 0) / r;
				mask = l > 0 ? 1 : 0;
				for (int j = 0; j < paramSize; j++)
				{
					d = l * value[j] * (1 - io);
					d0 = (1 - l)  * value[j] * io * mask;

					target[j][i] += d + d0;
				}
				xPos += hx;
			}
			yPos += hy;
		}
		zPos += hz;
	}
}

void StableFluid3d::TrilinearAdd(float x, float y, float z, float value, float * array)
{
	if (z > -maxZ && z < maxZ && y > -maxY && y < maxY && x > -maxX && x < maxX)
	{
		auto xPos = (x + maxX) / rangeX * resX;
		auto xFloor = (float)floor(xPos);
		auto dx = xPos - xFloor;
		auto mdx = 1 - dx;
		auto xId = (int)xFloor;
		auto xId1 = xId + 1;
		auto xBound = xId1 < resX;

		auto yPos = (y + maxY) / rangeY * resY;
		auto yFloor = (float)floor(yPos);
		auto dy = yPos - yFloor;
		auto mdy = 1 - dy;
		auto yId = (int)yFloor;
		auto yId1 = yId + 1;
		auto yBound = yId1 < resY;

		auto zPos = (z + maxZ) / rangeZ * resZ;
		auto zFloor = (float)floor(zPos);
		auto dz = zPos - zFloor;
		auto mdz = 1 - dz;
		auto zId = (int)zFloor;
		auto zId1 = zId + 1;
		auto zBound = zId1 < resZ;

		array[IX(zId, yId, xId)] += value * mdx * mdy * mdz;
		if (xBound)
			array[IX(zId, yId, xId1)] += value * dx * mdy * mdz;
		if (yBound)
		{
			array[IX(zId, yId1, xId)] += value * mdy * dy * mdz;
			if (xBound)
				array[IX(zId, yId1, xId1)] += value * dx * dy * mdz;
		}
		if (zBound)
		{
			array[IX(zId1, yId, xId)] += value * mdx * mdy * dz;
			if (xBound)
				array[IX(zId1, yId, xId1)] += value * dx * mdy * dz;
			if (yBound)
			{
				array[IX(zId1, yId1, xId)] += value * mdy * dy * dz;
				if (xBound)
					array[IX(zId1, yId1, xId1)] += value * dx * dy * dz;
			}
		}
	}
}

//--- private ----------------------

void StableFluid3d::AddSource(float * d, const float * s, const float dt)
{
	for (int l=0; l<length; l++)
		d[l] += s[l] * dt;
}

void StableFluid3d::CurlU()
{
	int id;
	float dw_dy, dv_dz;
	for (int z = 1; z <= nZ; z++)
	{
		for (int y = 1; y <= nY; y++)
		{
			for (int x = 1; x < avxX; x+=8)
			{
				id = IX(z, y, x);
				__m256 a = _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(_w + id + yf), _mm256_loadu_ps(_w + id - yf)), _mm256_set1_ps(0.5f));
				__m256 b = _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(_v + id + zf), _mm256_loadu_ps(_v + id - zf)), _mm256_set1_ps(0.5f));
				_mm256_storeu_ps(_cu + id, _mm256_sub_ps(a,b));
			}
			for (int x = avxX; x <= nX; x++)
			{
				id = IX(z, y, x);
				dw_dy = (_w[id + yf] - _w[id - yf]) * 0.5f;
				dv_dz = (_v[id + zf] - _v[id - zf]) * 0.5f;

				_cu[id] = dw_dy - dv_dz;
			}
		}
	}
}

void StableFluid3d::CurlV()
{
	int id;
	float du_dz, dw_dx;
	for (int z = 1; z <= nZ; z++)
	{
		for (int y = 1; y <= nY; y++)
		{
			for (int x = 1; x < avxX; x+=8)
			{
				id = IX(z, y, x);
				__m256 a = _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(_u + id + zf), _mm256_loadu_ps(_u + id - zf)), _mm256_set1_ps(0.5f));
				__m256 b = _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(_w + id +  1), _mm256_loadu_ps(_w + id -  1)), _mm256_set1_ps(0.5f));
				_mm256_storeu_ps(_cv + id, _mm256_sub_ps(a, b));
			}
			for (int x = avxX; x <= nX; x++)
			{
				id = IX(z, y, x);
				du_dz = (_u[id+zf] - _u[id-zf]) * 0.5f;
				dw_dx = (_w[id+1] - _w[id-1]) * 0.5f;

				_cv[id] = du_dz - dw_dx;
			}
		}
	}
}

void StableFluid3d::CurlW()
{
	int id;
	float dv_dx, du_dy;
	for (int z = 1; z <= nZ; z++)
	{
		for (int y = 1; y <= nY; y++)
		{
			for (int x = 1; x < avxX; x+=8)
			{
				id = IX(z, y, x);
				__m256 a = _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(_v + id +  1), _mm256_loadu_ps(_v + id -  1)), _mm256_set1_ps(0.5f));
				__m256 b = _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(_u + id + yf), _mm256_loadu_ps(_u + id - yf)), _mm256_set1_ps(0.5f));
				_mm256_storeu_ps(_cw + id, _mm256_sub_ps(a, b));
			}
			for (int x = avxX; x <= nX; x++)
			{
				id = IX(z, y, x);
				dv_dx = (_v[id+1] - _v[id-1]) * 0.5f;
				du_dy = (_u[id+yf] - _u[id-yf]) * 0.5f;

				_cw[id] = dv_dx - du_dy;
			}
		}
	}
}

void StableFluid3d::VorticityConfinement(float * u0, float * v0, float * w0)
{
	int id;
	float dw_dx, dw_dy, dw_dz;
	float l;
	float u, v, w;
	float res[4]{ 0,0,0,0 };

	for (int z = 2; z < nZ; z++)
	{
		for (int y = 2; y < nY; y++)
		{
			for (int x = 2; x < nX; x++)
			{
				id = IX(z, y, x);
				// Find derivative of the magnitude (n = del |w|)
				int xp = id + 1;
				int xm = id - 1;
				//dw_dx = (sqrt(_cu[xp] * _cu[xp] + _cv[xp] * _cv[xp] + _cw[xp] * _cw[xp]) -
				//		sqrt(_cu[xm] * _cu[xm] + _cv[xm] * _cv[xm] + _cw[xm] * _cw[xm])
				//		) * 0.5f;
				int yp = id + yf;
				int ym = id - yf;
				//dw_dy = (sqrt(_cu[yp] * _cu[yp] + _cv[yp] * _cv[yp] + _cw[yp] * _cw[yp]) -
				//		sqrt(_cu[ym] * _cu[ym] + _cv[ym] * _cv[ym] + _cw[ym] * _cw[ym])
				//		) * 0.5f;
				int zp = id + zf;
				int zm = id - zf;
				//dw_dz = (sqrt(_cu[zp] * _cu[zp] + _cv[zp] * _cv[zp] + _cw[zp] * _cw[zp]) -
				//		sqrt(_cu[zm] * _cu[zm] + _cv[zm] * _cv[zm] + _cw[zm] * _cw[zm])
				//		) * 0.5f;

				__m128 a = _mm_set_ps(_cu[xp], _cu[yp], _cu[zp], 1.0);
				a = _mm_mul_ps(a, a);
				__m128 b = _mm_set_ps(_cv[xp], _cv[yp], _cv[zp], 1.0);
				b = _mm_mul_ps(b, b);
				__m128 c = _mm_set_ps(_cw[xp], _cw[yp], _cw[zp], 1.0);
				c = _mm_mul_ps(c, c);
				
				__m128 sqrt1 = _mm_sqrt_ps(_mm_add_ps(a, _mm_add_ps(b, c)));

				a = _mm_set_ps(_cu[xm], _cu[ym], _cu[zm], 1.0);
				a = _mm_mul_ps(a, a);
				b = _mm_set_ps(_cv[xm], _cv[ym], _cv[zm], 1.0);
				b = _mm_mul_ps(b, b);
				c = _mm_set_ps(_cw[xm], _cw[ym], _cw[zm], 1.0);
				c = _mm_mul_ps(c, c);

				__m128 sqrt2 = _mm_sqrt_ps(_mm_add_ps(a, _mm_add_ps(b, c)));

				_mm_storeu_ps(res, _mm_mul_ps(_mm_set1_ps(0.5f), _mm_sub_ps(sqrt1, sqrt2)));
				dw_dx = res[3];
				dw_dy = res[2];
				dw_dz = res[1];
				// Calculate vector length. (|n|)
				// Add small factor to prevent divide by zeros.
				l = 1.0f / (sqrt(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001f);

				// N = ( n/|n| )
				dw_dx *= l;
				dw_dy *= l;
				dw_dz *= l;

				u = _cu[id];
				v = _cv[id];
				w = _cw[id];

				// N x w
				u0[id] = dw_dy * w - dw_dz * v;
				v0[id] = dw_dz * u - dw_dx * w;
				w0[id] = dw_dx * v - dw_dy * u;
			}
		}
	}
}

void StableFluid3d::SetBounds(const int b, float * d)
{
	for (int z = 1; z <= nZ; z++)
	{
		for (int y = 1; y <= nY; y++)
		{
			d[IX(z, y, 0)] = b == 1 ? -d[IX(z, y, 1)] : d[IX(z, y, 1)]; //left side
			d[IX(z, y, nX + 1)] = b == 1 ? -d[IX(z, y, nX)] : d[IX(z, y, nX)]; //right side
		}
	}
	__m256 m;
	__m256 n;
	for (int z = 1; z <= nZ; z++)
	{
		for (int x = 1; x < avxX; x+=8)
		{
			m = _mm256_loadu_ps(d + IX(z, 1, x));
			n = _mm256_loadu_ps(d + IX(z, nY, x));
			if (b == 2)
			{
				m = _mm256_mul_ps(m, _mm256_set1_ps(-1));
				n = _mm256_mul_ps(n, _mm256_set1_ps(-1));
			}
			_mm256_storeu_ps(d + IX(z, 0, x), m);
			_mm256_storeu_ps(d + IX(z, nY + 1, x), n);
		}
		for (int x = avxX; x <= nX; x++)
		{
			d[IX(z, 0, x)] = b == 2 ? -d[IX(z, 1, x)] : d[IX(z, 1, x)]; //bottom 
			d[IX(z, nY + 1, x)] = b == 2 ? -d[IX(z, nY, x)] : d[IX(z, nY, x)]; //top
		}
	}

	for (int y = 1; y <= nY; y++)
	{
		for (int x = 1; x < avxX; x+=8)
		{
			m = _mm256_loadu_ps(d + IX(1, y, x));
			n = _mm256_loadu_ps(d + IX(nZ, y, x));
			if (b == 3)
			{
				m = _mm256_mul_ps(m, _mm256_set1_ps(-1));
				n = _mm256_mul_ps(n, _mm256_set1_ps(-1));
			}
			_mm256_storeu_ps(d + IX(0, y, x), m);
			_mm256_storeu_ps(d + IX(nZ + 1, y, x), n);
		}
		for (int x = avxX; x <= nX; x++)
		{
			d[IX(0, y, x)] = b == 3 ? -d[IX(1, y, x)] : d[IX(1, y, x)]; //front 
			d[IX(nZ + 1, y, x)] = b == 3 ? -d[IX(nZ, y, x)] : d[IX(nZ, y, x)]; //back
		}
	}

	//edges
	for (int x = 1; x < avxX; x+=8)
	{
		_mm256_storeu_ps(d + x, _mm256_mul_ps(_mm256_set1_ps(0.5f), _mm256_add_ps(_mm256_loadu_ps(d + IX(0, 1, x)), _mm256_loadu_ps(d + IX(1, 0, x)))));
		_mm256_storeu_ps(d + IX(0, nY + 1, x), _mm256_mul_ps(_mm256_set1_ps(0.5f), _mm256_add_ps(_mm256_loadu_ps(d + IX(0, nY, x)), _mm256_loadu_ps(d + IX(1, nY + 1, x)))));
		_mm256_storeu_ps(d + IX(nZ + 1, 0, x), _mm256_mul_ps(_mm256_set1_ps(0.5f), _mm256_add_ps(_mm256_loadu_ps(d + IX(nZ + 1, 1, x)), _mm256_loadu_ps(d + IX(nZ, 0, x)))));
		_mm256_storeu_ps(d + IX(nZ + 1, nY + 1, x), _mm256_mul_ps(_mm256_set1_ps(0.5f), _mm256_add_ps(_mm256_loadu_ps(d + IX(nZ + 1, nY, x)), _mm256_loadu_ps(d + IX(nZ, nY + 1, x)))));
	}
	for (int x = avxX; x <= nX; x++)
	{
		d[IX(0, 0, x)] = 0.5f * (d[IX(0, 1, x)] + d[IX(1, 0, x)]); //bottom front
		d[IX(0, nY + 1, x)] = 0.5f * (d[IX(0, nY, x)] + d[IX(1, nY + 1, x)]); //top front
		d[IX(nZ + 1, 0, x)] = 0.5f * (d[IX(nZ + 1, 1, x)] + d[IX(nZ, 0, x)]); //bottom back
		d[IX(nZ + 1, nY + 1, x)] = 0.5f * (d[IX(nZ + 1, nY, x)] + d[IX(nZ, nY + 1, x)]); //top back
	}

	for (int y = 1; y <= nY; y++)
	{
		//left front
		d[IX(0, y, 0)] = 0.5f * (d[IX(0, y, 1)] + d[IX(1, y, 0)]);
		//right front
		d[IX(0, y, nX + 1)] = 0.5f * (d[IX(0, y, nX)] + d[IX(1, y, nX + 1)]);
		//left back
		d[IX(nZ + 1, y, 0)] = 0.5f * (d[IX(nZ + 1, y, 1)] + d[IX(nZ, y, 0)]);
		//right back
		d[IX(nZ + 1, y, nX + 1)] = 0.5f * (d[IX(nZ + 1, y, nX)] + d[IX(nZ, y, nX + 1)]);
	}

	for (int z = 1; z <= nZ; z++)
	{
		//left bottom 
		d[IX(z, 0, 0)] = 0.5f * (d[IX(z, 0, 1)] + d[IX(z, 1, 0)]);
		//right bottom
		d[IX(z, 0, nX + 1)] = 0.5f * (d[IX(z, 0, nX)] + d[IX(z, 1, nX + 1)]);
		//left top 
		d[IX(z, nY + 1, 0)] = 0.5f * (d[IX(z, nY + 1, 1)] + d[IX(z, nY, 0)]);
		//right top
		d[IX(z, nY + 1, nX + 1)] = 0.5f * (d[IX(z, nY + 1, nX)] + d[IX(z, nY, nX + 1)]);
	}

	//corners
	//left bottom front
	d[IX(0, 0, 0)] = 0.3333f * (d[IX(0, 0, 1)] + d[IX(0, 1, 0)] + d[IX(1, 0, 0)]);
	//right bottom front
	d[IX(0, 0, nX + 1)] = 0.3333f * (d[IX(0, 0, nX)] + d[IX(0, 1, nX + 1)] + d[IX(1, 0, nX + 1)]);
	//left top front
	d[IX(0, nY + 1, 0)] = 0.3333f * (d[IX(0, nY + 1, 1)] + d[IX(0, nY, 0)] + d[IX(1, nY + 1, 0)]);
	//right top front
	d[IX(0, nY + 1, nX + 1)] = 0.3333f * (d[IX(0, nY + 1, nX)] + d[IX(0, nY, nX + 1)] + d[IX(1, nY + 1, nX + 1)]);
	//left bottom back
	d[IX(nZ + 1, 0, 0)] = 0.3333f * (d[IX(nZ + 1, 0, 1)] + d[IX(nZ + 1, 1, 0)] + d[IX(nZ, 0, 0)]);
	//right bottom back
	d[IX(nZ + 1, 0, nX + 1)] = 0.3333f * (d[IX(nZ + 1, 0, nX)] + d[IX(nZ + 1, 1, nX + 1)] + d[IX(nZ, 0, nX + 1)]);
	//left top back
	d[IX(nZ + 1, nY + 1, 0)] = 0.3333f * (d[IX(nZ + 1, nY + 1, 1)] + d[IX(nZ + 1, nY, 0)] + d[IX(nZ, nY + 1, 0)]);
	//right top back
	d[IX(nZ + 1, nY + 1, nX + 1)] = 0.3333f * (d[IX(nZ + 1, nY + 1, nX)] + d[IX(nZ + 1, nY, nX + 1)] + d[IX(nZ, nY + 1, nX + 1)]);
}

void StableFluid3d::LinearSolve(const int b, float * d, const float * d0, const float a, const float c)
{
	int id;
	float ax = hxp*a;
	float ay = hyp*a;
	float az = hzp*a;
	float * swap;
	
	memcpy(_scratch, d, length * sizeof(float));
	for (int k = 0; k < SolverIterations; k++)
	{
		for (int z = 1; z <= nZ; z++)
		{
			for (int y = 1; y <= nY; y++)
			{
				for (int x = 1; x < avxX; x+=8)
				{
					id = IX(z, y, x);

					__m256 dx = _mm256_mul_ps(_mm256_add_ps(_mm256_loadu_ps(d + id -  1), _mm256_loadu_ps(d + id +  1)), _mm256_set1_ps(ax));
					__m256 dy = _mm256_mul_ps(_mm256_add_ps(_mm256_loadu_ps(d + id - yf), _mm256_loadu_ps(d + id + yf)), _mm256_set1_ps(ay));
					__m256 dz = _mm256_mul_ps(_mm256_add_ps(_mm256_loadu_ps(d + id - zf), _mm256_loadu_ps(d + id + zf)), _mm256_set1_ps(az));
					__m256 sum = _mm256_add_ps(_mm256_loadu_ps(d0 + id), _mm256_add_ps(dx, _mm256_add_ps(dy, dz)));
					_mm256_storeu_ps(_scratch + id, _mm256_div_ps(sum, _mm256_set1_ps(c)));
				}

				for (int x = avxX; x <= nX; x++)
				{
					id = IX(z, y, x);
					_scratch[id] = (d0[id] + ((d[id - 1] + d[id + 1]) *ax + (d[id - yf] + d[id + yf]) * ay + (d[id - zf] + d[id + zf]) * az)) / c;
				}
			}
		}
		SetBounds(b, _scratch);
		swap = _scratch;
		_scratch = d;
		d = swap;
	}
}

void StableFluid3d::Diffuse(const int b, float * d, const float * d0, const float diff, const float dt)
{
	float a = dt * diff;
	float c = 1 + dt * diff * 2 * (hxp + hyp + hzp);
	LinearSolve(b, d, d0, a, c);
}

void StableFluid3d::Advect(const int b, float * d, const float * d0, const float * u, const float * v, const float * w, const float dt)
{
	int id, id0, xId0, yId0, zId0, xId1, yId1, zId1;
	float q, r, s, dx, dy, dz, mdx, mdy, mdz;
	float dthx = dt / hx, dthy = dt / hy, dthz = dt / hz;

	float res[4]{ 0,0,0,0 };
	__m128 dth = _mm_set_ps(dthx, dthy, dthz, 0.1);
	__m128 mmin = _mm_set_ps(nX + 0.5f, nY + 0.5f, nZ + 0.5f, 1.0f);
	
	for (int z = 1; z <= nZ; z++)
	{
		for (int y = 1; y <= nY; y++)
		{
			for (int x = 1; x <= nX; x++)
			{
				id = IX(z, y, x);
				//q = x - u[id] * dthx;
				//q = fmax(0.5f, fmin(nX + 0.5f, q));
				//r = y - v[id] * dthy;
				//r = fmax(0.5f, fmin(nY + 0.5f, r));
				//s = z - w[id] * dthz;
				//s = fmax(0.5f, fmin(nZ + 0.5f, s));
				
				__m128 qrs = _mm_sub_ps(_mm_set_ps(x, y, z, 1.0f), _mm_mul_ps(_mm_set_ps(u[id], v[id], w[id], 1.0f), dth));
				_mm_storeu_ps(res, _mm_max_ps(_mm_set1_ps(0.5f), _mm_min_ps(qrs, mmin )));
				q = res[3]; //fuckr, little endian
				r = res[2];
				s = res[1];

				id0 = IX(s, r, q);
				xId0 = (int)q;
				xId1 = xId0 + 1;
				yId0 = (int)r;
				yId1 = yId0 + 1;
				zId0 = (int)s;
				zId1 = zId0 + 1;

				dx = q - xId0;
				mdx = 1 - dx;
				dy = r - yId0;
				mdy = 1 - dy;
				dz = s - zId0;
				mdz = 1 - dz;
				d[id] = mdx * mdy * mdz * d0[IX(zId0,yId0,xId0)] +
					dx * mdy * mdz * d0[IX(zId0, yId0, xId1)] +
					mdx * dy * mdz * d0[IX(zId0, yId1, xId0)] +
					dx * dy * mdz * d0[IX(zId0, yId1, xId1)] +
					mdx * mdy * dz * d0[IX(zId1, yId0, xId0)] +
					dx * mdy * dz * d0[IX(zId1, yId0, xId1)] +
					mdx * dy * dz * d0[IX(zId1, yId1, xId0)] +
					dx * dy * dz * d0[IX(zId1, yId1, xId1)];
				//d[id] = mdx * mdy * mdz * d0[id0] +
				//	dx * mdy * mdz * d0[id0 + 1] +
				//	mdx * dy * mdz * d0[id0 + yf] +
				//	dx * dy * mdz * d0[id0 + yf + 1] +
				//	mdx * mdy * dz * d0[id0 + zf] +
				//	dx * mdy * dz * d0[id0 + zf + 1] +
				//	mdx * dy * dz * d0[id0 + zf + yf] +
				//	dx * dy * dz * d0[id0 + zf + yf + 1];
			}
		}
	}
	SetBounds(b, d);
}

void StableFluid3d::Project(float * u, float * v, float * w, float * p, float * div)
{
	int id;
	for (int z = 1; z <= nZ; z++)
	{
		for (int y = 1; y <= nY; y++)
		{
			for (int x = 1; x < avxX; x+=8)
			{
				id = IX(z, y, x);
				__m256 ud = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(u + id +  1), _mm256_loadu_ps(u + id -  1)), _mm256_set1_ps(hx));
				__m256 vd = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(v + id + yf), _mm256_loadu_ps(v + id - yf)), _mm256_set1_ps(hy));
				__m256 wd = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(w + id + zf), _mm256_loadu_ps(w + id - zf)), _mm256_set1_ps(hz));
				__m256 divRes = _mm256_mul_ps(_mm256_set1_ps(-0.5f), _mm256_add_ps(ud, _mm256_add_ps(vd, wd)));
				_mm256_storeu_ps(div + id, divRes);
				_mm256_storeu_ps(p + id, _mm256_set1_ps(0.0f));
			}
			for (int x = avxX; x <= nX; x++)
			{
				id = IX(z, y, x);
				div[id] = -0.5f * ((u[id+ 1] - u[id- 1]) / hx +
					(v[id+yf] - v[id-yf]) / hy +
					(w[id+zf] - w[id-zf]) / hz);
				p[id] = 0;
			}
		}
	}
	SetBounds(0, div);
	SetBounds(0, p);

	float a = 1;
	float c = 2 * (hxp + hyp + hzp);
	LinearSolve(0, p, div, a, c);

	//subtract gradient
	for (int z = 1; z <= nZ; z++)
	{
		for (int y = 1; y <= nY; y++)
		{
			for (int x = 1; x < avxX; x+=8)
			{
				id = IX(z, y, x);
				__m256 ud   = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(p + id + 1), _mm256_loadu_ps(p + id - 1)), _mm256_set1_ps(hx));
				__m256 utmp = _mm256_loadu_ps(u + id);
				_mm256_storeu_ps(u + id, _mm256_sub_ps(utmp, _mm256_mul_ps(ud, _mm256_set1_ps(0.5f))));

				__m256 vd = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(p + id + yf), _mm256_loadu_ps(p + id - yf)), _mm256_set1_ps(hy));
				__m256 vtmp = _mm256_loadu_ps(v + id);
				_mm256_storeu_ps(v + id, _mm256_sub_ps(vtmp, _mm256_mul_ps(vd, _mm256_set1_ps(0.5f))));

				__m256 wd = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(p + id + zf), _mm256_loadu_ps(p + id - zf)), _mm256_set1_ps(hz));
				__m256 wtmp = _mm256_loadu_ps(w + id);
				_mm256_storeu_ps(w + id, _mm256_sub_ps(wtmp, _mm256_mul_ps(wd, _mm256_set1_ps(0.5f))));

				//__m128 ud = _mm_div_ps(_mm_sub_ps(_mm_loadu_ps(p + id + 1), _mm_loadu_ps(p + id - 1)), _mm_set1_ps(hx));
				//__m128 utmp = _mm_loadu_ps(u + id);
				//_mm_storeu_ps(u + id, _mm_sub_ps(utmp, _mm_mul_ps(ud, _mm_set1_ps(0.5f))));
			}
			for (int x = avxX; x <= nX; x++)
			{
				id = IX(z, y, x);
				u[id] -= 0.5f * (p[id + 1] - p[id - 1]) / hx;
				v[id] -= 0.5f * (p[id + yf] - p[id - yf]) / hy;
				w[id] -= 0.5f * (p[id + zf] - p[id - zf]) / hz;
			}
		}
	}
	/*for (int z = 1; z <= nZ; z++)
		for (int y = 1; y <= nY; y++)
			for (int x = 1; x <= nX; x++)
				u[IX(z, y, x)] -= 0.5f * (p[IX(z, y, x + 1)] - p[IX(z, y, x - 1)]) / hx;*/
	SetBounds(1, u);

	/*for (int z = 1; z <= nZ; z++)
		for (int y = 1; y <= nY; y++)
			for (int x = 1; x <= nX; x++)
				v[IX(z, y, x)] -= 0.5f * (p[IX(z, y + 1, x)] - p[IX(z, y - 1, x)]) / hy;*/
	SetBounds(2, v);
	
	/*for (int z = 1; z <= nZ; z++)
		for (int y = 1; y <= nY; y++)
			for (int x = 1; x <= nX; x++)
				w[IX(z, y, x)] -= 0.5f * (p[IX(z + 1, y, x)] - p[IX(z - 1, y, x)]) / hz;*/
	SetBounds(3, w);
}