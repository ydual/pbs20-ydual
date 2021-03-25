#pragma once

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include "constants.h"

using namespace Eigen;

// Particle Base class: Define several basic feature for Material Point Particle
class MPMParticle3D {
public:
	float v_p0;				// Initial Volume
	float m_p;				// Particle Mass

	Vector3f x_p;			// Particle Position
	Vector3f u_p;			// Particle Velocity
	Matrix3f B_p;			// Affine Momentum Matrix
	Matrix3f A_p;			// For computation purpose


	MPMParticle3D(const float in_v_p0, const float in_m_p, const Vector3f& in_x_p, const Vector3f& in_u_p, const Matrix3f& in_B_p) :
		v_p0(in_v_p0), m_p(in_m_p), x_p(in_x_p), u_p(in_u_p), B_p(in_B_p) {};
	// Member Functions
	virtual void ConstitutiveModel() {};										// Deformation gradient increment
	virtual void UpdateDeformation(const Matrix3f& T, float m_dt) {};				// Deformation gradient update
};

// Water class
class Water3D : public MPMParticle3D {
public:
	// Class Members
	//float A_p;
	float J_p;

	Water3D(const float in_v_p0, const float in_m_p, const Vector3f& in_x_p, const Vector3f& in_u_p, const Matrix3f& in_B_p) :
		MPMParticle3D(in_v_p0, in_m_p, in_x_p, in_u_p, in_B_p), J_p(1.0f) {
		A_p.setZero();
	};

	Water3D(MPMParticle3D& P) : MPMParticle3D(P.v_p0, P.m_p, P.x_p, P.u_p, P.B_p), J_p(1.0f) {
		A_p.setZero();
	};

	// Member Functions
	virtual void ConstitutiveModel();										// Deformation gradient increment
	virtual void UpdateDeformation(const Matrix3f& T, float m_dt);				// Deformation gradient update

};

// Sand class
class Sand3D : public MPMParticle3D {
public:
	/* Data */
	Matrix3f Fe, FeTr;											// (Trial) Elastic deformation
	Matrix3f Fp, FpTr;											// (Trial) Plastic deformation

	double q, alpha;												// Hardening paremeters

	double r;													// Color


	Sand3D(const float in_v_p0, const float in_m_p, const Vector3f& in_x_p, const Vector3f& in_u_p, const Matrix3f& in_B_p) : MPMParticle3D(in_v_p0, in_m_p, in_x_p, in_u_p, in_B_p) {
		A_p.setZero();

		Fe.setIdentity(); FeTr.setIdentity();
		Fp.setIdentity(); FpTr.setIdentity();

		q = 0.0;
		double phi = H0 + (H1 * q - H3) * exp(-H2 * q);
		alpha = sqrt(2 / 3.0) * 2 * sin(phi) / (3 - sin(phi));

		r = ((double)rand() / (RAND_MAX));
	}
	Sand3D(MPMParticle3D P) :MPMParticle3D(P.v_p0, P.m_p, P.x_p, P.u_p, P.B_p) {
		A_p.setZero();

		Fe.setIdentity(); FeTr.setIdentity();
		Fp.setIdentity(); FpTr.setIdentity();

		q = 0.0;
		double phi = H0 + (H1 * q - H3) * exp(-H2 * q);
		alpha = sqrt(2 / 3.0) * 2 * sin(phi) / (3 - sin(phi));

		r = ((double)rand() / (RAND_MAX));
	}

	/* Functions */
	virtual void ConstitutiveModel();								// Deformation gradient increment
	virtual void UpdateDeformation(const Matrix3f& T, float m_dt);				// Deformation gradient update
	void Plasticity();										// Update plastic dissipation
	void Projection(const Matrix3f& Eps, Matrix3f* T, float* dq);											// Return mapping algorithm

};

