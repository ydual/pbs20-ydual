#ifndef RIGIDOBJECTLS_H
#define RIGIDOBJECTLS_H

#include "RigidObject.h"

enum RigidObjectType{
    BALL, 
    CUBE
};

class RigidObjectLS: public RigidObject {
public:
    RigidObjectLS(): RigidObject() {}
    RigidObjectLS(const std::string &mesh_path, const ObjType t = ObjType::DYNAMIC):
        RigidObject(mesh_path, t) {}

    RigidObjectType ro_type;
    void setRigidObjectType(const RigidObjectType &t) {
        ro_type = t;
    }
    double getLevelSetValue(const Eigen::Vector3d& point) const;
    Eigen::Vector3d getLevelSetNormal(const Eigen::Vector3d& point) const;

    // temporal storage for forces and torques. In use if the force and torque are applied together after all the collisions within one time stamp.
    Eigen::Vector3d t_force; // momentum
    Eigen::Vector3d t_torque; // angular momentum
    void reset_for_particle_collision() {
        t_force = Eigen::Vector3d(0.0, 0.0, 0.0);
        t_torque = Eigen::Vector3d(0.0, 0.0, 0.0);
    }
    void add_t_force(const Eigen::Vector3d &v) {
        t_force += v;
    }
    void add_t_torque(const Eigen::Vector3d &v) {
        t_torque += v;
    }
    Eigen::Vector3d get_t_force() {
        return t_force;
    }
    Eigen::Vector3d get_t_torque() {
        return t_torque;
    }

private:
    Eigen::Vector3d getRelativePoint(const Eigen::Vector3d& point) const;
    // ball
    double getLevelSetValue_ball(const Eigen::Vector3d& point) const;
    Eigen::Vector3d getLevelSetNormal_ball(const Eigen::Vector3d& point) const;
    // cube
    double getLevelSetValue_cube(const Eigen::Vector3d& point) const;
    Eigen::Vector3d getLevelSetNormal_cube(const Eigen::Vector3d& point) const;
};

#endif
