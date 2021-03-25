#pragma once
#include <Eigen/Core>

using namespace Eigen;


/* ----- Changeable Parameters ----- */

// Grid size
const static int X_GRID = 200;
const static int Y_GRID = 100;
const static int Z_GRID = 200;

// // grid size for debugging
// const static int X_GRID = 100;
// const static int Y_GRID = 50;
// const static int Z_GRID = 100;

// // grid size for castle
// const static int X_GRID = 330;
// const static int Y_GRID = 300;
// const static int Z_GRID = 100;

// initilizatition source
const static int INIT_SOURCE = 0;

// APIC parameters
#define INTERPOLATION 1								// Kernel used in G2P and P2G: [1] Cubic - [2] Quadratic
const static float DT = 0.001;						    // Time-step

// Material
#define Material Sand									// [Water] - [Sand]
#define Material3D Sand3D								// [Sand3D] - [Water3D]

// Epsilong for float
const double eps = 1e-7;

/* ----- Grid ----- */
const static float H_INV = 1.0;

/* ----- Kernel Parameter ----- */
#if INTERPOLATION == 1
const static int CUB = 2;
const static Vector2f Translation_xp = Vector2f::Zero();
const static Vector3f Translation3D_xp = Vector3f::Zero();
static const int bni = -1;
static const float Dp_scal = 3.0;

#elif INTERPOLATION == 2
const static float CUB = 1.5;
const static Vector2f Translation_xp = Vector2f(0.5f, 0.5f);
static const int bni = 0;
static const float Dp_scal = 4.0;
#endif

/* ----- MATERIAL POINT PARAMETER ----- */

//#if Material == Water
//#define FRICTION false
//#else
//#define FRICTION true
//#endif
#define FRICTION true

const static Vector2f G = Vector2f(0.0f, -9.81);		// Gravity
const static Vector3f G_3D = Vector3f(0.0f, -9.81, 0.0f);		// Gravity3D (linpa)
const static float CFRI = 0.3;							// Friction coefficient	

/* Water */
const static float RHO_water = 1.0f;					// Density
const static float K_water = 50.0f;						// Bulk Modulus
const static int   GAMMA_water = 3;					// Penalize deviation form incompressibility

const static int DT_ROB = 30;

/* Sand */
static const float RHO_dry_sand = 2200.0f;
static const float E_dry_sand = 3.537e5;
static const float V_dry_sand = 0.3f;
static const float LAMBDA_dry_sand = E_dry_sand * V_dry_sand / (1.0f + V_dry_sand) / (1.0f - 2.0f * V_dry_sand);
static const float MU_dry_sand = E_dry_sand / (1.0f + V_dry_sand) / 2.0f;

static const float PI = 3.1415927f;
static const float H0 = 35 * PI / 180.0f;
static const float H1 = 9 * PI / 180.0f;
static const float H2 = 0.2f;
static const float H3 = 10 * PI / 180.0f;

/* Collision */
static const float collision_thres = 1.0;

// Kernel Functions
static float Bspline(float x)					// Cubic Bspline
{
	float W;
	x = fabs(x);

	if (x < 1)
		W = (x * x * x / 2.0f - x * x + 2 / 3.0f);

	else if (x < 2)
		W = (2 - x) * (2 - x) * (2 - x) / 6.0f;

	else
		W = 0;

	return W;
}

static float dBspline(float x)					// Cubic Bspline derivative
{
	float dW;
	float x_abs;
	x_abs = fabs(x);

	if (x_abs < 1)
		dW = 1.5f * x * x_abs - 2.0f * x;

	else if (x_abs < 2)
		dW = -x * x_abs / 2.0f + 2 * x - 2 * x / x_abs;

	else
		dW = 0;

	return dW;
}

static float getw_ip(const Vector2f& dist)		// 2D weight
{
	return Bspline(dist[0]) * Bspline(dist[1]);
}

static Vector2f getdw_ip(const Vector2f& dist)	// 2D weight gradient
{
	return Vector2f(
		dBspline(dist[0]) * Bspline(dist[1]),
		Bspline(dist[0]) * dBspline(dist[1]));
}

static float getw_ip(const Vector3f& dist)		// 3D weight
{
	return Bspline(dist[0]) * Bspline(dist[1]) * Bspline(dist[2]);
}

static Vector3f getdw_ip(const Vector3f& dist)	// 3D weight gradient
{
	return Vector3f(
		dBspline(dist[0]) * Bspline(dist[1]) * Bspline(dist[2]),
		Bspline(dist[0]) * dBspline(dist[1]) * Bspline(dist[2]),
		Bspline(dist[0]) * Bspline(dist[1]) * dBspline(dist[2]));
}
