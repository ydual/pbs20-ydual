#include "Grid3D.h"

bool Grid3d::NodeCollisionsRigidObjects(std::vector<RigidObjectLS>& m_objects, double epsilon, double scale, bool use_rotation) {
    u_i_co = u_i;
    // use_rotation = false;

    for (auto obj = m_objects.begin(); obj != m_objects.end(); ++obj) {
        Eigen::Vector3d disp = (x_i / scale).cast<double>() - obj->getPosition();
        if (obj->getLevelSetValue(disp) < 0) {
            Eigen::Vector3d n = obj->getLevelSetNormal(disp);
            Eigen::Vector3d vrel_vec;
            if (use_rotation)
                vrel_vec = (u_i / scale).cast<double>() - obj->getVelocity((x_i / scale).cast<double>());
            else
                vrel_vec = (u_i / scale).cast<double>() - obj->getLinearVelocity();
            double vrel = n.dot(vrel_vec);
            if (vrel > 0) {
                // bodies are moving apart
                return true; // one particle only collides with one object
            }
            Eigen::Vector3d numerator = - (1+epsilon) * vrel_vec;
            double t1 = 1.0 / m_i;
            double t2 = obj->getMassInv();
            double t4 = 0.0;
            if (use_rotation)
                t4 = n.dot((obj->getInertiaInvWorld() * (disp.cross(n))).cross(disp));
            Eigen::Vector3d j = numerator / (t1 + t2 + t4);
            Eigen::Vector3d force = j.dot(n) * n;

            u_i_co += scale * force.cast<float>() / m_i;
            obj->setLinearMomentum(obj->getLinearMomentum() - force);
            if (use_rotation)
                obj->setAngularMomentum(obj->getAngularMomentum() - disp.cross(force));
            return true; // one particle only collides with one object
        }
    }
    return false;
}

void Grid3d::NodeCollisions(std::vector<Boundary3D>& boundaries) {
	u_i_c = u_i_co;
	for (int i = 0; i < boundaries.size(); i++) {
		boundaries[i].Collision(x_i, u_i_c, collision_index, i);
	}
}

void Grid3d::NodeFrictions(std::vector<Boundary3D>& boundaries) {
	u_i_f = u_i_c;
	for (int i = 0; i < collision_index.size(); i++) {
		boundaries[collision_index[i]].Friction(u_i_f, u_i_c, u_i);
	}
};

void Grid3d::ResetGrid()
{
	m_i = 0.0f;
	u_i = Vector3f(0.0f, 0.0f, 0.0f);
	f_i = Vector3f(0.0f, 0.0f, 0.0f);
	collision_index.clear();
}
