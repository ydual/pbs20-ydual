#pragma once
#include "constants.h"
#include <vector>
using namespace Eigen;

// Boundary class handle collision
class Boundary {
public:
	Vector2f norm;
	int type;				// Sticky, Sperating or sliding

	std::vector<Vector2f> corners;

	Boundary(const Vector2f& in_norm, const int in_type, const std::vector<Vector2f>& in_corners) : norm(in_norm), type(in_type), corners(in_corners) {};

	void Collision(const Vector2f& node_coordinates, Vector2f& node_velocity, std::vector<int>& collision, const int b);// Apply collision with border object
	void Friction(Vector2f& Vi_fri, const Vector2f& Vi_col, const Vector2f& Vi);						// Apply friction with border object
		


};