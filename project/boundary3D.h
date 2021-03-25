#pragma once
#include "constants.h"
#include <vector>
#include <Eigen/Dense>
using namespace Eigen;

// Boundary class handle collision
class Boundary3D {
public:
	Vector3f norm;
	int type;				// Sticky, Sperating or sliding

	std::vector<Vector3f> corners;
    std::vector<double> min_vec, max_vec;
    int x_idx, y_idx; bool supported;

    Boundary3D(const Vector3f& in_norm, const int in_type, const std::vector<Vector3f>& in_corners) : norm(in_norm), type(in_type), corners(in_corners) {
        min_vec.resize(3);
        std::fill(min_vec.begin(), min_vec.end(), std::numeric_limits<double>::max());
        max_vec.resize(3);
        std::fill(max_vec.begin(), max_vec.end(), std::numeric_limits<double>::min());
        for (auto corner: in_corners) {
            if (corner[0] > max_vec[0])
                max_vec[0] = corner[0];
            if (corner[0] < min_vec[0])
                min_vec[0] = corner[0];
            if (corner[1] > max_vec[1])
                max_vec[1] = corner[1];
            if (corner[1] < min_vec[1])
                min_vec[1] = corner[1];
            if (corner[2] > max_vec[2])
                max_vec[2] = corner[2];
            if (corner[2] < min_vec[2])
                min_vec[2] = corner[2];
        }
        x_idx = 0; y_idx = 1; supported = true;
	    if (abs(abs(norm[0]) - 1) < eps) {
	    	x_idx = 1;
	    	y_idx = 2;
	    }
	    else if (abs(abs(norm[1]) - 1) < eps) {
	    	x_idx = 0;
	    	y_idx = 2;
	    }
	    else if (abs(abs(norm[2]) - 1) < eps) {
	    	x_idx = 0;
	    	y_idx = 1;
	    }
        else
            supported = false;
    };

	bool Collision(const Vector3f& node_coordinates, Vector3f& node_velocity, std::vector<int>& collision, const int &b);// Apply collision with border object
	void Friction(Vector3f& Vi_fri, const Vector3f& Vi_col, const Vector3f& Vi);						// Apply friction with border object
	bool isInside(const Vector3f& projectedPos); // ensure the projected position is inside the boundary formed by corners
};

