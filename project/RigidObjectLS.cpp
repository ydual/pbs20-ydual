#include "RigidObjectLS.h"
#include <iostream>
#include <cmath>

Eigen::Vector3d RigidObjectLS::getRelativePoint(const Eigen::Vector3d &point) const {
    Eigen::Matrix3d R = getRotationMatrix();
    Eigen::Vector3d p_object = getPosition();
    Eigen::Vector3d p = R.transpose() * (point - p_object);
    return p;
}

double RigidObjectLS::getLevelSetValue(const Eigen::Vector3d& point) const {
    switch(ro_type) {
        case BALL: return getLevelSetValue_ball(point);
        case CUBE: return getLevelSetValue_cube(point);
    }
    throw std::runtime_error("Not implemented");
}

Eigen::Vector3d RigidObjectLS::getLevelSetNormal(const Eigen::Vector3d& point) const{
    switch(ro_type) {
        case BALL: return getLevelSetNormal_ball(point);
        case CUBE: return getLevelSetNormal_cube(point);
    }
    throw std::runtime_error("Not implemented");
}

double RigidObjectLS::getLevelSetValue_ball(const Eigen::Vector3d& point) const{
    auto p = getRelativePoint(point);
    double scale = getScale();
    double radius = scale * scale * 19.7371;
    return p.norm() - radius;
}

Eigen::Vector3d RigidObjectLS::getLevelSetNormal_ball(const Eigen::Vector3d& point) const{
    auto p = getRelativePoint(point);
    Eigen::Matrix3d R = getRotationMatrix();
    return R * (p / p.norm());
}

double RigidObjectLS::getLevelSetValue_cube(const Eigen::Vector3d& point) const{
    auto p = getRelativePoint(point);
    double scale = getScale();
    double radius = scale * scale * 1.0;
    double x, y, z;
    x = p[0]; y = p[1]; z = p[2];
    double abs_x, abs_y, abs_z;
    abs_x = std::abs(x); abs_y = std::abs(y); abs_z = std::abs(z);
    return std::max(std::max(abs_x, abs_y), abs_z) - radius;
}

Eigen::Vector3d RigidObjectLS::getLevelSetNormal_cube(const Eigen::Vector3d& point) const{
    auto p = getRelativePoint(point);
    Eigen::Matrix3d R = getRotationMatrix();
    double x, y, z;
    x = p[0]; y = p[1]; z = p[2];
    double abs_x, abs_y, abs_z;
    abs_x = std::abs(x); abs_y = std::abs(y); abs_z = std::abs(z);
    if (abs_x > abs_y && abs_x > abs_z)
        return R * Eigen::Vector3d(x/abs_x, 0.0, 0.0);
    else if (abs_y > abs_x && abs_y > abs_z)
        return R * Eigen::Vector3d(0.0, y/abs_y, 0.0);
    else if (abs_z > abs_x && abs_z > abs_y)
        return R * Eigen::Vector3d(0.0, 0.0, z/abs_z);
    else if (abs_x == abs_y && abs_x > abs_z)
        return R * Eigen::Vector3d(x/abs_x, y/abs_y, 0.0) / sqrt(2);
    else if (abs_x == abs_z && abs_x > abs_y)
        return R * Eigen::Vector3d(x/abs_x, 0.0, z/abs_z) / sqrt(2);
    else if (abs_y == abs_z && abs_y > abs_x)
        return R * Eigen::Vector3d(0.0, y/abs_y, z/abs_z) / sqrt(2);
    else // (abs_x == abs_y && abs_x == abs_z)
        return R * Eigen::Vector3d(x/abs_x, y/abs_y, z/abs_z) / sqrt(3);
}

