#pragma once

/// \brief basic fluid simulation in 2d space
class StableFluid2d
{
private:
	int length; 	///< total number of elements in simulation grid
	int resX, nX; 	///< number of elements in x dimension
	int resY, nY; 	///< number of elements in y dimension

	float rangeX, hx, hxp, maxX; ///< simulation space size in x dimensions
	float rangeY, hy, hyp, maxY; ///< simulation space size in y dimensions

	float * _u; ///< main velocity x component array
	float * _v; ///< main velocity y component array
	float * _u0; ///< secondary velocity x component array
	float * _v0; ///< secondary velocity y component array

	float * _dens; ///< main density array
	float * _dens0; ///< secondary density array
			  
	float * _curl; ///< curl array for vorticity confinement

	float * _scratch; ///< scratch array for the solver

	/// SSE / AVX stuff
	int yf;
	int sseX;
	int avxX;

	/// \brief helper function translating 2d indexing to linear 
	inline int IX(int y, int x) { return x + y*resX; };

	/// \brief add values to a simulation component
	/// \param d target component to be added to
	/// \param s values to be added
	/// \param dt simulation timestep
	void AddSource(float * d, const float * s, const float & dt);

	/// \brief avoids dissipation of velocity details
	/// \param u0 x velocity component target to add to
	/// \param v0 y velocity component target to add to
	/// \param curl scratch array for curl value
	/// \param u source x velocity component 
	/// \param v source y velocity component 
	void VorticityConfinement(float * u0, float * v0, float * curl, const float * u, const float * v);

	/// \brief enforce boundary conditions on component array
	/// \param b indicates usage: 0 for any, [1..n] spatial dimension
	/// \param d component to be adjusted
	void SetBounds(const int b, float * d);

	/// \brief Jacobi solver for stable diffusion
	/// \param b indicates usage: 0 for any, [1..n] spatial dimension
	/// \param d target component
	/// \param d0 source component
	void LinearSolve(const int b, float * d, const float * d0, const float a, const float c);

	/// \brief diffusion step: interaction of neighbouring values
	/// \param b indicates usage: 0 for any, [1..n] spatial dimension
	/// \param d target component
	/// \param d0 source component
	/// \param diff diffusion rate
	/// \param dt simulation timestep
	void Diffuse(const int b, float * d, const float * d0, const float diff, const float & dt);

	/// \brief advection step: actual moving of values by a linear backtrace
	/// \param b indicates usage: 0 for any, [1..n] spatial dimension
	/// \param d target component
	/// \param d0 source component
	/// \param u source x velocity component
	/// \param u source y velocity component
	/// \param dt simulation timestep
	void Advect(const int b, float * d, const float * d0, const float * u, const float * v, const float & dt);

	/// \brief mass conserving routine for velocities by calculating and subtracting the gradient
	/// \param u target x velocity component
	/// \param v target y velocity component
	/// \param p pressure (i guess) scratch array
	/// \param div divergence (i guess) scratch array
	void Project(float* u, float* v, float* p, float* div);

public:
	/// \brief constructor
	/// \param countX resolution in x dimension
	/// \param countY resolution in y dimension
	/// \param width size in x dimension
	/// \param height size in y dimension
	StableFluid2d(int countX, int countY, float width, float height);

	/// \brief destructor
	~StableFluid2d();

	/// \brief solver iteration count
	int SolverIterations = 20;
			 
	int Count() { return length; };			///< returns the total number of elements in simulation grid
	int CountX() { return resX; };			///< returns the total number of elements in x dimension
	int CountY() { return resY; };			///< returns the total number of elements in y dimension
	float RangeX() { return rangeX; };		///< returns simulation space size in x dimension
	float RangeY() { return rangeY; };		///< returns simulation space size in y dimension
	float HX() { return hx; };				///< returns the grid spacing in x dimension
	float HY() { return hy; };				///< returns the grid spacing in y dimension

	float* U() { return _u; };					///< returns velocity x component values
	float* V() { return _v; };					///< returns velocity y component values
	float* Density() { return _dens; };			///< returns density values
	float* NewU() { return _u0; };				///< returns velocity x component values to be added
	float* NewV() { return _v0; };				///< returns velocity y component values to be added
	float* NewDensity() { return _dens0; };		///< returns density component values to be added

	/// \brief resets the arrays holding the values to be added to 0
	void ResetInputs();

	/// \brief updates the simulation
	/// \param visc viscosity amount
	/// \param diff diffusion amount
	/// \param simulation timestep
	void Update(const float visc, const float diff, const float dt);

	/// \brief multiplies all velocities by a certain rate
	void MultiplyVelocity(const float rate);

	/// \brief multiplies all densities by a certain rate
	void MultiplyDensity(const float rate);

	/// \brief adds a value to an array around a given position
	/// \param cx center x of insert position
	/// \param cy center y of insert position
	/// \param sx scale factor of the radius in x dimension
	/// \param sy scale factor of the radius in y dimension
	/// \param r radius around the center
	/// \param io lerps the value along the radius between [value..0] to [0..value]
	/// \param paramSize count of values and their corrensponding target arrays
	/// \param value array of values to be added
	/// \param target array of target arrays to which the values should be added
	void AddRadial(const float cx, const float cy, const float sx, const float sy, const float r, const float io, const int paramSize, const float * value, float ** target);
};