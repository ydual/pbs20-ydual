#include "Grid2D.h"

void Grid2d::NodeCollisions(std::vector<Boundary>& boundaries) {
	u_i_c = u_i;
	for (int i = 0; i < boundaries.size(); i++) {
		boundaries[i].Collision(x_i, u_i_c, collision_index, i);
	}
}

void Grid2d::NodeFrictions(std::vector<Boundary>& boundaries) {
	u_i_f = u_i_c;
	for (int i = 0; i < collision_index.size(); i++) {
		boundaries[collision_index[i]].Friction(u_i_f, u_i_c, u_i);
	}
};

void Grid2d::ResetGrid()
{
	m_i = 0.0f;
	u_i = Vector2f(0.0f, 0.0f);
	f_i = Vector2f(0.0f, 0.0f);
	collision_index.clear();
}
