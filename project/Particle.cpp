#include "Particle.h"

void Water::ConstitutiveModel() {
	float d_Jp = -K_water * (1.0f / pow(J_p, GAMMA_water) - 1.0f);
	A_p = d_Jp * v_p0 * J_p;
}

void Water::UpdateDeformation(const Matrix2f& T, float m_dt)
{
	J_p = (1 + m_dt * T.trace()) * J_p;
}

void Sand::ConstitutiveModel() {

	JacobiSVD<Matrix2f> svd(Fe, ComputeFullU | ComputeFullV);
	Matrix2f Eps = svd.singularValues().asDiagonal();

	// We have to pass through array() because log() not available on Matrices
	Matrix2f logEps = Eps.diagonal().array().log().matrix().asDiagonal();
	Matrix2f U = svd.matrixU();
	Matrix2f V = svd.matrixV();

	Matrix2f dFe =
		2 * MU_dry_sand * Eps.inverse() * logEps + LAMBDA_dry_sand * logEps.trace() * Eps.inverse();

	A_p = v_p0 * U * dFe * V.transpose() * Fe.transpose();
}

void Sand::UpdateDeformation(const Matrix2f& T, float m_dt)
{
	FeTr = (Matrix2f::Identity() + m_dt * T) * Fe;
	FpTr = Fp;

	Sand::Plasticity();
}

void Sand::Plasticity()
{
	JacobiSVD<Matrix2f> svd(FeTr, ComputeFullU | ComputeFullV);		// SVD decomposition
	Matrix2f Eps = svd.singularValues().asDiagonal();
	Matrix2f U = svd.matrixU();
	Matrix2f V = svd.matrixV();

	Matrix2f T; float dq;
	Sand::Projection(Eps, &T, &dq);

	Fe = U * T * V.transpose();
	Fp = V * T.inverse() * Eps * V.transpose() * FpTr;

	// Hardening
	q += dq;
	float phi = H0 + (H1 * q - H3) * exp(-H2 * q);
	alpha = (float)(sqrt(2.0 / 3.0) * (2.0 * sin(phi)) / (3.0f - sin(phi)));
}

void Sand::Projection(const Matrix2f& Eps, Matrix2f* T, float* dq)
{
	Matrix2f e, e_c;

	e = Eps.diagonal().array().log().matrix().asDiagonal();
	e_c = e - e.trace() / 2.0f * Matrix2f::Identity();

	if (e_c.norm() < 1e-8 || e.trace() > 0) {
		T->setIdentity();
		*dq = e.norm();
		return;										// Projection to the tip of the cone (case 2)
	}

	float dg = e_c.norm()
		+ (LAMBDA_dry_sand + MU_dry_sand) / MU_dry_sand * e.sum() * alpha;
	if (dg <= 0) {
		*T = Eps;
		*dq = 0;
		return;										// No projection (case 1)
	}

	Matrix2f Hm = e - dg * e_c / e_c.norm();

	*T = Hm.diagonal().array().exp().matrix().asDiagonal();
	*dq = dg;
	return;											// Projection onto the yield surface (case 3)
}

