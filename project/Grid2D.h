#pragma once

#include <iostream>
#include <cmath>

#include "constants.h"
#include "boundary.h"

using namespace Eigen;

class Grid2d {
public:
	float m_i;			// Mass of a grid node

	Vector2f x_i;
	Vector2f u_i;
	Vector2f u_i_c;			// velocity after collision
	Vector2f u_i_f;			// velocity after friction
	Vector2f f_i;

	std::vector<int> collision_index;			// List of boundary objects that collides the grid

	Grid2d(const Vector2f& in_x_i) : x_i(in_x_i) {
		m_i = 0.0f;
		u_i.setZero();
		u_i_c.setZero();
		u_i_f.setZero();

		f_i.setZero();
	};

	void NodeCollisions(std::vector<Boundary>& boundaries);								// Apply collision to all borders
	void NodeFrictions(std::vector<Boundary>& boundaries);								// Apply friction if collision

	void ResetGrid();

};
