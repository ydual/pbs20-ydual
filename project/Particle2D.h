#pragma once

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include "constants.h"

using namespace Eigen;

// Particle Base class: Define several basic feature for Material Point Particle
class MPMParticle {
public:
	float v_p0;				// Initial Volume
	float m_p;				// Particle Mass

	Vector2f x_p;			// Particle Position
	Vector2f u_p;			// Particle Velocity
	Matrix2f B_p;			// Affine Momentum Matrix

	MPMParticle(const float in_v_p0, const float in_m_p, const Vector2f& in_x_p, const Vector2f& in_u_p, const Matrix2f& in_B_p) :
		v_p0(in_v_p0), m_p(in_m_p), x_p(in_x_p), u_p(in_u_p), B_p(in_B_p) {};
};

// Water class
class Water : public MPMParticle {
public:
	// Class Members
	float A_p;
	float J_p;

	Water(const float in_v_p0, const float in_m_p, const Vector2f& in_x_p, const Vector2f& in_u_p, const Matrix2f& in_B_p) :
		MPMParticle(in_v_p0, in_m_p, in_x_p, in_u_p, in_B_p), A_p(0.0f), J_p(1.0f) {};

	Water(const MPMParticle& P) : MPMParticle(P.v_p0, P.m_p, P.x_p, P.u_p, P.B_p), A_p(0.0f), J_p(1.0f) {};

	// Member Functions
	void ConstitutiveModel();										// Deformation gradient increment
	void UpdateDeformation(const Matrix2f& T, float m_dt);				// Deformation gradient update

};

// Sand class
class Sand : public MPMParticle {
public:
	/* Data */
	Matrix2f A_p;												// For computation purpose

	Matrix2f Fe, FeTr;											// (Trial) Elastic deformation
	Matrix2f Fp, FpTr;											// (Trial) Plastic deformation

	double q, alpha;												// Hardening paremeters
	double r;													// Color

	Sand(const float in_v_p0, const float in_m_p, const Vector2f& in_x_p, const Vector2f& in_u_p, const Matrix2f& in_B_p) : MPMParticle(in_v_p0, in_m_p, in_x_p, in_u_p, in_B_p) {
		A_p.setZero();

		Fe.setIdentity(); FeTr.setIdentity();
		Fp.setIdentity(); FpTr.setIdentity();

		q = 0.0;
		double phi = H0 + (H1 * q - H3) * exp(-H2 * q);
		alpha = sqrt(2 / 3.0) * 2 * sin(phi) / (3 - sin(phi));

		r = ((double)rand() / (RAND_MAX));
	}
	Sand(const MPMParticle &P) :MPMParticle(P.v_p0, P.m_p, P.x_p, P.u_p, P.B_p) {
		A_p.setZero();

		Fe.setIdentity(); FeTr.setIdentity();
		Fp.setIdentity(); FpTr.setIdentity();

		q = 0.0;
		double phi = H0 + (H1 * q - H3) * exp(-H2 * q);
		alpha = sqrt(2 / 3.0) * 2 * sin(phi) / (3 - sin(phi));

		r = ((double)rand() / (RAND_MAX));
	}

	/* Functions */
	void ConstitutiveModel();								// Deformation gradient increment
	void UpdateDeformation(const Matrix2f& T, float m_dt);				// Deformation gradient update
	void Plasticity();										// Update plastic dissipation
	void Projection(const Matrix2f& Eps, Matrix2f* T, float* dq);											// Return mapping algorithm

};
