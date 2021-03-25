#include "Particle3D.h"

void Water3D::ConstitutiveModel() {
	float d_Jp = -K_water * (1.0f / pow(J_p, GAMMA_water) - 1.0f);
	A_p = d_Jp * v_p0 * J_p * Matrix3f::Identity();
}

void Water3D::UpdateDeformation(const Matrix3f& T, float m_dt)
{
	J_p = (1 + m_dt * T.trace()) * J_p;
}

// 3D version sand
void Sand3D::ConstitutiveModel() {

	JacobiSVD<Matrix3f> svd(Fe, ComputeFullU | ComputeFullV);
	Matrix3f Eps = svd.singularValues().asDiagonal();

	// We have to pass through array() because log() not available on Matrices
	Matrix3f logEps = Eps.diagonal().array().log().matrix().asDiagonal();
	Matrix3f U = svd.matrixU();
	Matrix3f V = svd.matrixV();

	Matrix3f dFe =
		2 * MU_dry_sand * Eps.inverse() * logEps + LAMBDA_dry_sand * logEps.trace() * Eps.inverse();

	A_p = v_p0 * U * dFe * V.transpose() * Fe.transpose();
}

void Sand3D::UpdateDeformation(const Matrix3f& T, float m_dt)
{
	FeTr = (Matrix3f::Identity() + m_dt * T) * Fe;
	FpTr = Fp;

	Sand3D::Plasticity();
}

void Sand3D::Plasticity()
{
	JacobiSVD<Matrix3f> svd(FeTr, ComputeFullU | ComputeFullV);		// SVD decomposition
	Matrix3f Eps = svd.singularValues().asDiagonal();
	Matrix3f U = svd.matrixU();
	Matrix3f V = svd.matrixV();

	Matrix3f T; float dq;
	Sand3D::Projection(Eps, &T, &dq);

	Fe = U * T * V.transpose();
	Fp = V * T.inverse() * Eps * V.transpose() * FpTr;

	// hardening
	q += dq;
	float phi = H0 + (H1 * q - H3) * exp(-H2 * q);
	alpha = (float)(sqrt(2.0 / 3.0) * (2.0 * sin(phi)) / (3.0f - sin(phi)));
}

void Sand3D::Projection(const Matrix3f& Eps, Matrix3f* T, float* dq)
{
	Matrix3f e, e_c;

	e = Eps.diagonal().array().log().matrix().asDiagonal();
	e_c = e - e.trace() / 3.0f * Matrix3f::Identity(); // 3D version

	if (e_c.norm() < 1e-8 || e.trace() > 0) {
		T->setIdentity();
		*dq = e.norm();
		return;										// Projection to the tip of the cone (case 2)
	}

	float dg = e_c.norm()
		+ (3.0 * LAMBDA_dry_sand + 2.0 * MU_dry_sand) / (2.0 * MU_dry_sand) * e.sum() * alpha; // dimenmsion changed

	if (dg <= 0) {
		*T = Eps;
		*dq = 0;
		return;										// No projection (case 1)
	}

	Matrix3f Hm = e - dg * e_c / e_c.norm();

	*T = Hm.diagonal().array().exp().matrix().asDiagonal();
	*dq = dg;
	return;											// Projection onto the yield surface (case 3)
}


