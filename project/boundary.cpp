#include "boundary.h"

void Boundary::Collision(const Vector2f& node_coordinates, Vector2f& node_velocity, std::vector<int>& collision, const int b)
{
	// Collide when the point is inside the boundary
	float distance = norm.dot(node_coordinates - corners[0]);


	//Sticky collisions enforce that a point remains fixed to a
	//particular reference point on the collision object.
	if (type == 1 && distance < 0) {
		node_velocity = Vector2f(0.0f, 0.0f);
	}

	// Seperating: If a node is already inside a collision body,
	// then it should not penetrate any deeper.
	// If a node is originally outside the object, then it should remain outside

	Vector2f next_position = node_coordinates + DT * node_velocity;
	double next_distance = norm.dot(next_position - corners[0]);
	double dist_c = next_distance - std::min(distance, 0.0f);

	if (type == 2 && dist_c < 0) {
		node_velocity -= dist_c * norm / DT;
		collision.push_back(b);
	}

	// Slipping:we do not want to allow separation for existing collisions, 
	// but sliding along the surface is permitted.

	if (type == 3 && distance < 0) {
		node_velocity -= dist_c * norm / DT;
		collision.push_back(b);
	}
}

void Boundary::Friction(Vector2f& Vi_fri, const Vector2f& Vi_col, const Vector2f& Vi) {
	// Compute Tangential Velocity
	Vector2f V_t = Vi_col - norm * (norm.dot(Vi_fri));

	if (V_t.norm() > eps)
	{

		Vector2f t = V_t / V_t.norm();

		// Coulomb's friction 
		Vi_fri -= std::min(V_t.norm(), CFRI * (Vi_col - Vi).norm()) * t;
	}
}