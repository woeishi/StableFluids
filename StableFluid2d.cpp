#include "stdafx.h"
#include "StableFluid2d.h"

#include <algorithm>

StableFluid2d::StableFluid2d(int countX, int countY, float width, float height)
{
	length = countX * countY;
	resX = countX;
	nX = countX - 2;
	resY = countY;
	nY = countY - 2;

	rangeX = width;
	hx = rangeX / countX;
	hxp = pow(hx, -2.0f);
	maxX = (rangeX - hx) / 2.0f;

	rangeY = height;
	hy = rangeY / countY;
	hyp = pow(hy, -2.0f);
	maxY = (rangeY - hy) / 2.0f;

	_u = new float[length]();
	_v = new float[length]();
	_u0 = new float[length]();
	_v0 = new float[length]();

	_dens = new float[length]();
	_dens0 = new float[length]();

	_curl = new float[length]();

	_scratch = new float[length]();

	// SSE/AVX
	yf = resX;
	auto dr = std::div(nX, 4); //1-nX inclusive
	sseX = dr.quot * 4 + 1; //1+ because inclusive

	dr = std::div(nX, 8);
	avxX = dr.quot * 8 + 1;
}

StableFluid2d::~StableFluid2d()
{
	delete[] _u, _u0;
	delete[] _v, _v0;
	delete[] _dens, _dens0;

	delete[] _curl;
	delete[] _scratch;
}

void StableFluid2d::ResetInputs()
{
	memset(_dens0, 0, length * sizeof(float));
	memset(_u0, 0, length * sizeof(float));
	memset(_v0, 0, length * sizeof(float));
}

void StableFluid2d::Update(const float visc, const float diff, const float dt)
{
	//velocity
	AddSource(_u, _u0, dt);
	AddSource(_v, _v0, dt);

	VorticityConfinement(_u0, _v0, _curl, _u, _v);
	AddSource(_u, _u0, dt);
	AddSource(_v, _v0, dt);

	Diffuse(1, _u0, _u, visc, dt);
	Diffuse(2, _v0, _v, visc, dt);
	Project(_u0, _v0, _u, _v);

	Advect(1, _u, _u0, _u0, _v0, dt);
	Advect(2, _v, _v0, _u0, _v0, dt);
	Project(_u, _v, _u0, _v0);

	//density
	AddSource(_dens, _dens0, dt);
	Diffuse(0, _dens0, _dens, diff, dt);
	Advect(0, _dens, _dens0, _u, _v, dt);

	ResetInputs();
}

void StableFluid2d::MultiplyVelocity(const float rate)
{
	for (int i = 0; i < length; i++)
	{
		_u[i] *= rate;
		_v[i] *= rate;
	}
}

void StableFluid2d::MultiplyDensity(const float rate)
{
	for (int i = 0; i < length; i++)
		_dens[i] *= rate;
}

void StableFluid2d::AddRadial(const float cx, const float cy, const float sx, const float sy, const float r, const float io, const int paramSize, const float * value, float ** target)
{
	//int paramSize = std::min(value.size(), target.size());
	auto dr = div(resX, 8);
	int max8 = dr.quot * 8;
	float off[8]{ 0,1,2,3,4,5,6,7 };
	int i = 0;
	float px, py, l, mask, d, d0;

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
			__m256 _l = _mm256_sqrt_ps(_mm256_add_ps(_px, _py));
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
			l = sqrt(px*px + py*py);
			l = fmax(r - l, 0) / r;
			mask = l > 0 ? 1 : 0;

			for (int j = 0; j < paramSize; j++)
			{
				d = l * value[j] * (1 - io);
				d0 = (1 - l) * value[j] * io * mask;

				target[j][i] += d + d0;
			}
			xPos += hx;
		}
		yPos += hy;
	}
}

//--- private ----------------------

void StableFluid2d::AddSource(float * d, const float * s, const float & dt)
{
	for (int l = 0; l < length; l++) //automatic vectorization is getting this loop
		d[l] += s[l] * dt;
}

void StableFluid2d::VorticityConfinement(float * u0, float * v0, float * curl, const float * u, const float * v)
{
	float du_dy, dv_dx;
	float dw_dx, dw_dy;
	float length;
	float c;
	int id;

	// Calculate magnitude of curl(u,v) for each cell. (|w|)
	for (int y = 1; y <= nY; y++)
	{
		for (int x = 1; x < avxX; x += 8)
		{
			id = IX(y, x);
			__m256 a = _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(v + id + 1), _mm256_loadu_ps(v + id - 1)), _mm256_set1_ps(0.5f));
			__m256 b = _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(u + id + yf), _mm256_loadu_ps(u + id - yf)), _mm256_set1_ps(0.5f));
			_mm256_storeu_ps(_curl + id, _mm256_sub_ps(a, b));
		}
		for (int x = avxX; x <= nX; x++)
		{
			id = IX(y, x);
			dv_dx = (v[id + 1] - v[id - 1]) * 0.5f;
			du_dy = (u[id + resX] - u[id - resX]) * 0.5f;
			_curl[id] = abs(du_dy - dv_dx);
		}
	}

	for (int y = 2; y < nY; y++)
	{
		for (int x = 2; x < nX; x++)
		{
			id = IX(y, x);
			// Find derivative of the magnitude (n = del |w|)
			dw_dx = (_curl[id + 1] - _curl[id - 1]) * 0.5f;
			dw_dy = (_curl[id + resX] - _curl[id - resX]) * 0.5f;

			// Calculate vector length. (|n|)
			// Add small factor to prevent divide by zeros.
			length = sqrt(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001f;

			// N = ( n/|n| )
			dw_dx /= length;
			dw_dy /= length;

			du_dy = (_u[id + resX] - _u[id - resX]) * 0.5f;
			dv_dx = (_v[id + 1] - _v[id - 1]) * 0.5f;
			c = du_dy - dv_dx;

			// N x w
			u0[id] = dw_dy * -c;
			v0[id] = dw_dx * c;
		}
	}
}

void StableFluid2d::SetBounds(const int b, float* d)
{
	int id = resX;
	int right = resX - 1;
	int top = length - resX;

	for (int y = resX; y < top; y += resX)
	{
		d[y] = b == 1 ? -d[y + 1] : d[y + 1]; //left side
		d[right + y] = b == 1 ? -d[right + y - 1] : d[right + y - 1]; //right side
	}

	__m256 m;
	__m256 n;
	for (int x = 1; x < avxX; x += 8)
	{
		m = _mm256_loadu_ps(d + x + resY);
		n = _mm256_loadu_ps(d + IX(nY, x));
		if (b == 2)
		{
			m = _mm256_mul_ps(m, _mm256_set1_ps(-1));
			n = _mm256_mul_ps(n, _mm256_set1_ps(-1));
		}
		_mm256_storeu_ps(d + x, m);
		_mm256_storeu_ps(d + IX(nY + 1, x), n);
	}
	for (int x = avxX; x <= nX; x++)
	{
		d[x] = b == 2 ? -d[x + resY] : d[x + resY]; //bottom 
		d[x + top] = b == 2 ? -d[x + top - resX] : d[x + top - resX]; //top
		//d[IX(nY + 1, x)] = b == 2 ? -d[IX(nY, x)] : d[IX(nY, x)]; //top
	}
	
	//left bottom
	d[0] = 0.5f * (d[1] + d[resY]);
	//left top
	d[IX(nY + 1, 0)] = 0.5f * (d[IX(nY + 1, 1)] + d[IX(nY, 0)]);
	//right bottom
	d[nX + 1] = 0.5f * (d[nX] + d[IX(1, nX + 1)]);
	//right top
	d[IX(nY + 1, nX + 1)] = 0.5f * (d[IX(nY + 1, nX)] + d[IX(nY, nX + 1)]);
}

//Gauss-Seidel doesn't work for vectorization due to backwards reference
//void StableFluid2d::LinearSolve(int b, float * d, float * d0, float a, float c)
//{
//	int id;
//	float ax = _hxp * a;
//	float ay = _hyp * a;
//	auto dr = div(nX, 4);
//	int xMax = dr.quot * 4 + 1;
//	for (int k = 0; k < 20; k++)
//	{
//		for (int y = 1; y <= nY; y++)
//		{
//			for (int x = 1; x <= nX; x++)
//			{
//				//d[IX(y, x)] = (d0[IX(y, x)] + (d[IX(y, x - 1)] * _hxp*a + d[IX(y, x + 1)] * _hxp*a + d[IX(y - 1, x)] * _hyp*a + d[IX(y + 1, x)] * _hyp*a)) / c;
//				id = IX(y, x);
//				d[id] = (d0[id] + (d[id-1] * _hxp*a + d[id+1] * _hxp*a + d[id-resX] * _hyp*a + d[id+resX] * _hyp*a)) / c;
//			}
//		}
//		SetBounds(b, d);
//	}
//}

void StableFluid2d::LinearSolve(const int b, float * d, const float * d0, const float a, const float c)
{
	int id;
	float ax = hxp * a;
	float ay = hyp * a;
	
	memcpy(_scratch, d, length*sizeof(float));
	for (int k = 0; k < SolverIterations; k++)
	{
		for (int y = 1; y <= nY; y++)
		{
			for (int x = 1; x < avxX; x += 8)
			{
				id = IX(y, x);

				__m256 dx = _mm256_mul_ps(_mm256_add_ps(_mm256_loadu_ps(d + id - 1), _mm256_loadu_ps(d + id + 1)), _mm256_set1_ps(ax));
				__m256 dy = _mm256_mul_ps(_mm256_add_ps(_mm256_loadu_ps(d + id - resX), _mm256_loadu_ps(d + id + resX)), _mm256_set1_ps(ay));
				__m256 sum = _mm256_add_ps(_mm256_loadu_ps(d0 + id), _mm256_add_ps(dx, dy));
				__m256 factor = _mm256_div_ps(sum, _mm256_set1_ps(c));
				_mm256_storeu_ps(_scratch + id, factor);
			}

			for (int x = avxX; x <= nX; x++)
			{
				id = IX(y, x);
				_scratch[id] = (d0[id] + ((d[id - 1] + d[id + 1]) * ax  + (d[id - resX] + d[id + resX]) * ay)) / c;
			}
		}
		SetBounds(b, _scratch);
		std::swap(_scratch, d);
	}
}

void StableFluid2d::Diffuse(const int b, float * d, const float * d0, const float diff, const float & dt)
{
	float a = dt * diff;
	float c = 1 + dt * diff * 2 * (hxp + hyp);
	LinearSolve(b, d, d0, a, c);
}

void StableFluid2d::Advect(const int b, float * d, const float * d0, const float * u, const float * v, const float & dt)
{
	int id, xId0, yId0, xId1, yId1;
	float q, r, dx, dy, mdx, mdy;
	for (int y = 1; y <= nY; y++)
	{
		for (int x = 1; x <= nX; x++)
		{
			id = IX(y, x);
			q = x - u[id] * dt / hx;
			q = fmax(0.5f, fmin(nX + 0.5f, q));
			r = y - v[id] * dt / hy;
			r = fmax(0.5f, fmin(nY + 0.5f, r));

			xId0 = (int)q;
			xId1 = xId0 + 1;
			yId0 = (int)r;
			yId1 = yId0 + 1;

			dx = q - xId0;
			mdx = 1 - dx;
			dy = r - yId0;
			mdy = 1 - dy;
			d[id] = mdx * (mdy * d0[IX(yId0, xId0)] + dy * d0[IX(yId1, xId0)]) +
				dx * (mdy * d0[IX(yId0, xId1)] + dy * d0[IX(yId1, xId1)]);
		}
	}
	SetBounds(b, d);
}

void StableFluid2d::Project(float * u, float * v, float * p, float * div)
{
	//compute gradient == height field
	//solving Poisson equation
	int id;
	for (int y = 1; y <= nY; y++)
	{
		for (int x = 1; x < avxX; x+=8)
		{
			id = IX(y, x);
			__m256 ud = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(u + id + 1), _mm256_loadu_ps(u + id - 1)), _mm256_set1_ps(hx));
			__m256 vd = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(v + id + resX), _mm256_loadu_ps(v + id - resX)), _mm256_set1_ps(hy));
			__m256 divRes = _mm256_mul_ps(_mm256_set1_ps(-0.5f),_mm256_add_ps(ud,vd));
			_mm256_storeu_ps(div + id, divRes);
			_mm256_storeu_ps(p + id, _mm256_set1_ps(0.0f));
		}
		for (int x = avxX; x <= nX; x++)
		{
			id = IX(y, x);
			div[id] = -0.5f * ((u[id + 1] - u[id - 1]) / hx + (v[id + resX] - v[id - resX]) / hy);
			p[id] = 0;
		}
	}
	SetBounds(0, div); SetBounds(0, p);

	float a = 1;
	float c = 2 * (hxp + hyp);
	LinearSolve(0, p, div, a, c);

	//subtract gradient
	for (int y = 1; y <= nY; y++)
	{
		for (int x = 1; x < avxX; x+=8)
		{
			id = IX(y, x);
			__m256 ud = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(p + id + 1), _mm256_loadu_ps(p + id - 1)), _mm256_set1_ps(hx));
			__m256 utmp = _mm256_loadu_ps(u + id);
			_mm256_storeu_ps(u + id, _mm256_sub_ps(utmp,_mm256_mul_ps(ud,_mm256_set1_ps(0.5f))));

			__m256 vd = _mm256_div_ps(_mm256_sub_ps(_mm256_loadu_ps(p + id + resX), _mm256_loadu_ps(p + id - resX)), _mm256_set1_ps(hy));
			__m256 vtmp = _mm256_loadu_ps(v + id);
			_mm256_storeu_ps(v + id, _mm256_sub_ps(vtmp, _mm256_mul_ps(vd, _mm256_set1_ps(0.5f))));
		}
		for (int x = avxX; x <= nX; x++)
		{
			id = IX(y, x);
			u[id] -= 0.5f * (p[id + 1] - p[id - 1]) / hx;
			v[id] -= 0.5f * (p[id + resX] - p[id - resX]) / hy;
		}
	}
	SetBounds(1, u); SetBounds(2, v);
}
