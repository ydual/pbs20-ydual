#include "MPMSim3D.h"
#include <igl/Timer.h>

bool MPMSim3D::advance() {
    // advance
    advance_mpm(); // in the mpm phase, rigid bodies may also be updated via collisions.
    advance_rigid_bodies();
    
    // advance time
    m_time += m_dt;
    m_step++;
    return false;
}

bool MPMSim3D::advance_mpm() {
    // advance mpm 3d
    ENABLE_TIMER = false;
igl::Timer timer;
if (ENABLE_TIMER) {
    timer.start();
    std::cout<<"----------------------------------"<<std::endl;
    std::cout<<"Hey! Timer is starting!!!"<<std::endl;
}
    assert(int(valid_grid_indexes.size()) == 0);
    if (SCENEMODE == CITY || SCENEMODE == SANDPILE) {
        if (m_step % 5 == 0 && m_step < 5000) AddParticles();
    } 
    //else if (SCENEMODE == WATERMIXTURE) {
    //    if (m_step % 5 == 0 && m_step < 5000) AddParticles_water();
    //}
    
    P2G();
if (ENABLE_TIMER)
    std::cout<<"p2g time = \t\t"<<timer.getElapsedTime()<<" (n = "<<particles.size()<<")"<<std::endl;

    GetValidGrids();
if (ENABLE_TIMER)
    std::cout<<"get valid grid time = \t"<<timer.getElapsedTime()<<" (n = "<<grids.size()<<")"<<std::endl;

    UpdateGrids();
if (ENABLE_TIMER)
    std::cout<<"grid update time = \t"<<timer.getElapsedTime()<<" (n = "<<valid_grid_indexes.size()<<")"<<std::endl;

    G2P();
if (ENABLE_TIMER)
    std::cout<<"g2p time = \t\t"<<timer.getElapsedTime()<<" (n = "<<particles.size()<<")"<<std::endl;

    UpdateParticles();
if (ENABLE_TIMER)
    std::cout<<"particle update time = \t"<<timer.getElapsedTime()<<" (n = "<<particles.size()<<")"<<std::endl;

    ResetGrid();
if (ENABLE_TIMER)
    std::cout<<"grid reset time = \t"<<timer.getElapsedTime()<<" (n = "<<valid_grid_indexes.size()<<")"<<std::endl;

    valid_grid_indexes.clear();
    return false;
}

bool MPMSim3D::advance_rigid_bodies() {
    if (SCENEMODE == CITY) {
        // check for rigid body collision with the city
        for (auto o = m_objects_ls.begin(); o != m_objects_ls.end(); ++o) {
            if (o->ro_type != CUBE) continue;  // only CUBE case implemented
            for (int i = 0; i < 8; i++) {
                double cubScale = o->getScale() * o->getScale();
                Vector3d pos = Vector3d(cubScale - 2 * (i / 4) * cubScale, cubScale - 2 * ((i % 4) > 1) * cubScale, cubScale - 2 * (i % 2) * cubScale);
                Eigen::Vector3d disp = o->getRotationMatrix() * pos;
                Vector3d wordlPos = o->getRotationMatrix() * pos + o->getPosition();

                int x, y, z;
                double scale = 50.0;
                double epsilon = 0.9;
                x = wordlPos[0] * scale; y = wordlPos[1] * scale; z = wordlPos[2] * scale;
                y = std::max(0, y);
                if (x < 0 || z < 0) continue;
                else if (x >= X_GRID || y >= Y_GRID || z >= Z_GRID) continue;
                int idx = (X_GRID + 1) * (Y_GRID + 1) * z + (X_GRID + 1) * y + x;
                if (valueField[idx] < 0) {
                    Eigen::Vector3d n = normalField[idx].cast<double>();
                    Eigen::Vector3d vrel_vec;

                    vrel_vec = o->getVelocity(wordlPos);
                    double vrel = n.dot(vrel_vec);
                    if (vrel > 0) {
                        break; // one particle only collides with one object
                    }

                    Eigen::Vector3d numerator = -(1 + epsilon) * vrel_vec;
                    double t1 = 0; // consider the city motionless
                    double t2 = o->getMassInv();
                    double t4 = n.dot((o->getInertiaInvWorld() * (disp.cross(n))).cross(disp));
                    Eigen::Vector3d j = numerator / (t1 + t2 + t4);
                    Eigen::Vector3d force = j.dot(n) * n;

                    o->setLinearMomentum(o->getLinearMomentum() + force);
                    o->setAngularMomentum(o->getAngularMomentum() + (disp).cross(force));
                    break; // one particle only collides with one object
                }
            }
        }
    }

    // For now we ignore the collisions between rigid bodies, since re-implementing that broad-phase and narrow-phase detection is really annoying. 
    for (auto o = m_objects_ls.begin(); o != m_objects_ls.end(); ++o) {
        auto v = o->getLinearVelocity();
        auto p = o->getPosition();
        // if the object is on the floor, assume the force balance
        if (o->getPosition()[1] - o->getScale() * o->getScale() - CUB / double(50) < eps) {
           o->setLinearVelocity(v); // ignore gravity for the convenience
        }
        else {
            o->setLinearVelocity(v + DT * m_gravity);
        }
        o->setPosition(p + DT*o->getLinearVelocity());

        // update orientation
        Eigen::Vector3d w = o->getAngularVelocity();
        Eigen::Quaterniond wq;
        wq.w() = 0;
        wq.vec() = w;
        Eigen::Quaterniond q = o->getRotation();
        Eigen::Quaterniond dq = wq * q;
        Eigen::Quaterniond new_q;
        new_q.w() = q.w() + 0.5 * DT * dq.w();
        new_q.vec() = q.vec() + 0.5 * DT * dq.vec();
        o->setRotation(new_q.normalized());
    }
    return false;
}

void MPMSim3D::GetValidGrids() {
    const int num_threads = 64;
    int num_grids = grids.size();
    int num_grids_per_thread = std::ceil(float(num_grids) / num_threads);
    std::vector<std::vector<int>> valid_collections(num_threads);
#pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) {
        int idx_start = i * num_grids_per_thread;
        int idx_end = std::min((i + 1) * num_grids_per_thread, num_grids);
        for (int j = idx_start; j < idx_end; ++j)
            if (grids[j].m_i > 0)
                valid_collections[i].push_back(j);
    }
    for (int i = 0; i < num_threads; ++i) {
        for (auto j: valid_collections[i])
            valid_grid_indexes.push_back(j);
    }
}

void MPMSim3D::P2G() {
    int num_particles = particles.size();
#pragma omp parallel for 
    for (int i = 0; i < num_particles; i++) {
        // Pre-update Ap (in particle loop)
        particles[i]->ConstitutiveModel();
        // Index of bottom-left node closest to the particle
        int node_base = (X_GRID + 1) * (Y_GRID + 1) * max(min(static_cast<int>(particles[i]->x_p[2] - Translation3D_xp[2]), Z_GRID - CUB), CUB)
            + (X_GRID + 1) * max(min(static_cast<int>(particles[i]->x_p[1] - Translation3D_xp[1]), Y_GRID - CUB), CUB) + max(min(static_cast<int>(particles[i]->x_p[0] - Translation3D_xp[0]), X_GRID - CUB), CUB);


        // Loop over all the close nodes (depend on interpolation through bni)
        for (int z = bni; z < 3; z++) {
            for (int y = bni; y < 3; y++) {
                for (int x = bni; x < 3; x++) {
                    // Index of the node
                    int node_id = node_base + x + (X_GRID + 1) * y + (X_GRID + 1)* (Y_GRID + 1) * z;

                    // Distance and weight
                    Vector3f dist = particles[i]->x_p - grids[node_id].x_i;

                    float w_ip = getw_ip(dist);
                    Vector3f dw_ip = getdw_ip(dist);

                    // Pre-compute node mass, node velocity and pre-update force increment (APIC)
                    float in_m_i = w_ip * particles[i]->m_p;
                    Vector3f in_u_i = w_ip * particles[i]->m_p * (particles[i]->u_p + Dp_scal * H_INV * H_INV * H_INV * particles[i]->B_p * (-dist)); 
                    Vector3f in_f_i = particles[i]->A_p * dw_ip;

                    // Udpate mass, velocity and force
                    grids[node_id].m_i += in_m_i;
                    grids[node_id].u_i += in_u_i;
                    grids[node_id].f_i += in_f_i;
                }
            }
        }
    }
};

void MPMSim3D::UpdateGrids() {
igl::Timer timer;
if (ENABLE_TIMER) {
    timer.start();
}
    int num_valid_grids = valid_grid_indexes.size();
#pragma omp parallel for schedule (dynamic)
    for (int idx = 0; idx < num_valid_grids; idx++) {
        int i = valid_grid_indexes[idx];
        // Finish updating velocity, force, and apply updated force
        // velocity here is actually the momentum at first
        grids[i].u_i /= grids[i].m_i;
        grids[i].f_i = DT * (-grids[i].f_i / grids[i].m_i + G_3D);
        grids[i].u_i += grids[i].f_i;
    }
if (ENABLE_TIMER)
    std::cout<<"[SUB] grid update time = \t\t"<<timer.getElapsedTime()<<" (n = "<<particles.size()<<")"<<std::endl;

    int counter = 0;
    for (int idx = 0; idx < num_valid_grids; idx++) {
        int i = valid_grid_indexes[idx];
        if (grids[i].NodeCollisionsRigidObjects(m_objects_ls))
            counter++; // collide with rigid objects. Note that the status of rigid objects also change here.
    }
if (ENABLE_TIMER)
    std::cout<<"[SUB] collision with objects = \t\t"<<timer.getElapsedTime()<<" (n = "<<particles.size()<<")"<<std::endl;
    //std::cout<<"DEBUG: number of collisions: "<<counter<<std::endl;
    
#pragma omp parallel for schedule (dynamic)
    for (int idx = 0; idx < num_valid_grids; idx++) {
        int i = valid_grid_indexes[idx];
        grids[i].NodeCollisions(boundaries);
    }
if (ENABLE_TIMER)
    std::cout<<"[SUB] node collision with boundaries = \t"<<timer.getElapsedTime()<<" (n = "<<particles.size()<<")"<<std::endl;

    if (SCENEMODE == CITY)
        cityCollision(); 

#pragma omp parallel for schedule (dynamic)
    for (int idx = 0; idx < num_valid_grids; idx++) {
        int i = valid_grid_indexes[idx];
#if FRICTION
        grids[i].NodeFrictions(boundaries);
#else
        grids[i].u_i_f = grids[i].u_i_c;
#endif
    }

if (ENABLE_TIMER)
    std::cout<<"[SUB] node friction with boundaries = \t"<<timer.getElapsedTime()<<" (n = "<<particles.size()<<")"<<std::endl;

#if FRICTION
    if (SCENEMODE == CITY)
        cityFriction();
#endif

};

void MPMSim3D::G2P() {
    int num_particles = particles.size();
#pragma omp parallel for 
    for (int i = 0; i < num_particles; i++) {
        // Index of bottom-left node closest to the particle
        int node_base = (X_GRID + 1) * (Y_GRID + 1) * max(min(static_cast<int>(particles[i]->x_p[2] - Translation3D_xp[2]), Z_GRID - CUB), CUB)
            + (X_GRID + 1) * max(min(static_cast<int>(particles[i]->x_p[1] - Translation3D_xp[1]), Y_GRID - CUB), CUB) + max(min(static_cast<int>(particles[i]->x_p[0] - Translation3D_xp[0]), X_GRID - CUB), CUB);

        // Set velocity and velocity field to 0 for sum update
        particles[i]->u_p.setZero();
        particles[i]->B_p.setZero();

        // Loop over all the close nodes (depend on interpolation through bni)
        for (int z = bni; z < 3; z++) {
            for (int y = bni; y < 3; y++) {
                for (int x = bni; x < 3; x++) {
                    // Index of the node
                    int node_id = node_base + x + (X_GRID + 1) * y + (X_GRID + 1) * (Y_GRID + 1) * z; 

                    // Distance and weight
                    Vector3f dist = particles[i]->x_p - grids[node_id].x_i;
                    float w_ip = getw_ip(dist);

                    // Update velocity and velocity field (APIC)
                    particles[i]->u_p += w_ip * grids[node_id].u_i_f;
                    particles[i]->B_p += w_ip * (grids[node_id].u_i_f * (-dist).transpose());
                }
            }
        }

    }
};

void MPMSim3D::UpdateParticles() {
    int num_particles = particles.size();
#pragma omp parallel for 
    for (int i = 0; i < num_particles; i++) {
        // Index of bottom-left node closest to the particle
        int node_base = (X_GRID + 1) * (Y_GRID + 1) * max(min(static_cast<int>(particles[i]->x_p[2] - Translation3D_xp[2]), Z_GRID - CUB), CUB)
            + (X_GRID + 1) * max(min(static_cast<int>(particles[i]->x_p[1] - Translation3D_xp[1]), Y_GRID - CUB), CUB) + max(min(static_cast<int>(particles[i]->x_p[0] - Translation3D_xp[0]), X_GRID - CUB), CUB);

        // Save position to compute nodes-particle distances and update position in one loop
        Vector3f x_p_buff = particles[i]->x_p;

        //  T ~ nodal deformation
        Matrix3f T; T.setZero();
        Vector3f x_p_temp;
        x_p_temp.setZero();

        // Loop over all the close nodes (depend on interpolation through bni))
        for (int z = bni; z < 3; z++) {
            for (int y = bni; y < 3; y++) {
                for (int x = bni; x < 3; x++) {
                    // Index of the node
                    int node_id = node_base + x + (X_GRID + 1) * y + (X_GRID + 1) * (Y_GRID + 1) * z; 


                    // Distance and weight
                    Vector3f dist = x_p_buff - grids[node_id].x_i;
                    float w_ip = getw_ip(dist);
                    Vector3f dw_ip = getdw_ip(dist);

                    // Update position and nodal deformation
                    x_p_temp += w_ip * (grids[node_id].x_i + DT * grids[node_id].u_i_c);
                    T += grids[node_id].u_i_c * dw_ip.transpose();

                }
            }
        }
        particles[i]->x_p = x_p_temp;
        //if (particles[i]->x_p.norm() < eps) std::cout << "real 0" << std::endl;

        // Update particle deformation gradient (elasticity, plasticity etc...)
        particles[i]->UpdateDeformation(T, DT);
    }
};

void MPMSim3D::ResetGrid() {
    int num_valid_grids = valid_grid_indexes.size();
#pragma omp parallel for schedule (dynamic)
    for (int idx = 0; idx < num_valid_grids; idx++) {
        int i = valid_grid_indexes[idx];
        grids[i].ResetGrid();
    }
};

void MPMSim3D::cityCollision() {
    int num_valid_grids = valid_grid_indexes.size();

    int counter = 0;
#pragma omp parallel for schedule (dynamic)
    for (int idx = 0; idx < num_valid_grids; idx++) {
        int i = valid_grid_indexes[idx];
        int x = grids[i].x_i[0];
        int y = grids[i].x_i[1];
        int z = grids[i].x_i[2];
        int nodeIdx = (X_GRID + 1) * (Y_GRID + 1) * z + (X_GRID + 1) * y + x;
        if (valueField[nodeIdx] < 0) { // collision occurs
            counter++;
            float dist_c = normalField[nodeIdx].dot(DT * grids[i].u_i_c);

            grids[i].u_i_c -= dist_c * normalField[nodeIdx] / DT;
        }
    }
}

void MPMSim3D::cityFriction() {
    int num_valid_grids = valid_grid_indexes.size();

    int counter = 0;
#pragma omp parallel for schedule (dynamic)
    for (int idx = 0; idx < num_valid_grids; idx++) {
        int i = valid_grid_indexes[idx];
        int x = grids[i].x_i[0];
        int y = grids[i].x_i[1];
        int z = grids[i].x_i[2];
        int nodeIdx = (X_GRID + 1) * (Y_GRID + 1) * z + (X_GRID + 1) * y + x;
        if (valueField[nodeIdx] < 0) { // collision occurs
            Vector3f V_t = grids[i].u_i_c - normalField[nodeIdx] * (normalField[nodeIdx].dot(grids[i].u_i_f));
            if (V_t.norm() > eps) {
                Vector3f t = V_t / V_t.norm();
                grids[i].u_i_f -= std::min(V_t.norm(), CFRI * (grids[i].u_i_c - grids[i].u_i).norm()) * t;

            }
        }
    }
}

#undef ENABLE_MPM_TIMER

