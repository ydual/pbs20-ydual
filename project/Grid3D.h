#pragma once

#include <iostream>
#include <cmath>

#include "constants.h"
#include "boundary3D.h"
#include "RigidObjectLS.h"

using namespace Eigen;

class Grid3d {
public:
	float m_i;			// Mass of a grid node

	Vector3f x_i;		// position
	Vector3f u_i;		// velocity
	Vector3f u_i_co;		// velocity after collision with rigid objects
	Vector3f u_i_c;			// velocity after collision
	Vector3f u_i_f;			// velocity after friction
	Vector3f f_i;

	std::vector<int> collision_index;			// List of boundary objects that collides the grid

	Grid3d(const Vector3f& in_x_i) : x_i(in_x_i) {
		m_i = 0.0f;
		u_i.setZero();
		u_i_c.setZero();
		u_i_f.setZero();

		f_i.setZero();
	};

    bool NodeCollisionsRigidObjects(std::vector<RigidObjectLS>& m_objects, double epsilon=0.5, double scale=25.0, bool use_rotation=true);  // collision with rigid objects. Note that the status of the objects are also affected. 
	void NodeCollisions(std::vector<Boundary3D>& boundaries);								// Apply collision to all borders
	void NodeFrictions(std::vector<Boundary3D>& boundaries);								// Apply friction if collision

	void ResetGrid();

};
