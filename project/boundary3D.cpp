#include "boundary3D.h"
#include <iostream>

bool Boundary3D::Collision(const Vector3f& node_coordinates, Vector3f& node_velocity, std::vector<int>& collision, const int &b)
{
	// Collide when the point is inside the boundary
    float distance, dist_c;
    distance = norm.dot(node_coordinates - corners[0]); // particle "below" the surface
	double next_distance = norm.dot(node_coordinates + DT * node_velocity - corners[0]);

	if ((type == 1 || type == 3) && distance >= 0) {
		return false;
	}
    else if (type == 2) {
	    double next_distance = norm.dot(node_coordinates + DT * node_velocity - corners[0]);
	    double dist_c = next_distance - std::min(distance, 0.0f);
		if (dist_c >= 0)
			return false;
    }

	// check that the projected point is still inside the boundary
	Vector3f projectedPos = node_coordinates - distance * norm;
    if (!isInside(projectedPos)) 
        return false;

	//Sticky collisions enforce that a point remains fixed to a
	//particular reference point on the collision object.
	if (type == 1) {
		node_velocity = Vector3f(0.0f, 0.0f, 0.0f);
	}

	dist_c = next_distance - std::min(distance, 0.0f);
	// Seperating: If a node is already inside a collision body,
	// then it should not penetrate any deeper.
	// If a node is originally outside the object, then it should remain outside
	if (type == 2) {
		node_velocity -= dist_c * norm / DT;
		collision.push_back(b);
	}

	// Slipping:we do not want to allow separation for existing collisions, 
	// but sliding along the surface is permitted.
	if (type == 3) {
		node_velocity -= dist_c * norm / DT;
		collision.push_back(b);
	}
    return true;
}

void Boundary3D::Friction(Vector3f& Vi_fri, const Vector3f& Vi_col, const Vector3f& Vi) {
	// Compute Tangential Velocity
	Vector3f V_t = Vi_col - norm * (norm.dot(Vi_fri));

	if (V_t.norm() > eps)
	{

		Vector3f t = V_t / V_t.norm();
		// Coulomb's friction 
		Vi_fri -= std::min(V_t.norm(), CFRI * (Vi_col - Vi).norm()) * t;
	}
}

// works only when normal is in x/y/z-direction
bool Boundary3D::isInside(const Vector3f& projectedPos) {
	if (!supported) {
		std::cout << "not supported type of face" << std::endl;
		return false;
	}
    if (projectedPos[x_idx] > max_vec[x_idx] || projectedPos[x_idx] < min_vec[y_idx])
        return false;
    if (projectedPos[y_idx] > max_vec[y_idx] || projectedPos[y_idx] < min_vec[y_idx])
        return false;
	if (corners.size() < 2) return false;

	// we do not consider points on the boundary
	Matrix3f mat;
	mat << 1, corners[0][x_idx], corners[0][y_idx],
		1, corners[1][x_idx], corners[1][y_idx],
		1, projectedPos[x_idx], projectedPos[y_idx];
	bool sign = mat.determinant() < 0;

	// check that the point is on the same side as other corners
	for (size_t i = 1; i < corners.size(); i++) {
		Matrix3f new_mat;
		new_mat << 1, corners[i][x_idx], corners[i][y_idx],
			1, corners[(i + 1) % corners.size()][x_idx], corners[(i + 1) % corners.size()][y_idx],
			1, projectedPos[x_idx], projectedPos[y_idx];
        bool new_sign = new_mat.determinant() < 0;
		if (sign != new_sign) {
			return false;
		}
	}
	return true;
}

