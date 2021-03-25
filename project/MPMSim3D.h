#pragma once
#include "Simulation.h"
#include "Particle3D.h"
#include "constants.h"
#include "Grid3D.h"
#include "boundary.h"
#define POISSON_PROGRESS_INDICATOR 1
#include "PoissonGenerator.h"
#include "RigidObjectLS.h"

#include <igl/readOBJ.h>
#include <igl/readSTL.h>

enum MPMSimSceneType {
    CITY, // default
    CASTLE,
    TESTSCENE,
    SANDPILE,
    WATER
};

using namespace std;

// Different test scenes for our MPM solver
// Modified from the cannon ball exercise
class MPMSim3D : public Simulation {
public:
    MPMSimSceneType SCENEMODE = CITY;
    MPMSim3D(const int &scene_mode) : Simulation() {
        switch (scene_mode) {
            case 1: SCENEMODE = CITY; break;
            case 2: SCENEMODE = CASTLE; break;
            case 3: SCENEMODE = TESTSCENE; break;
            case 4: SCENEMODE = SANDPILE; break;
            case 5: SCENEMODE = WATER; break;
            default: SCENEMODE = CITY; break;
        }
        init();
    }

    virtual void init() override {
        m_objects_ls.clear();
        m_dt = 0.01f;
        m_mass = 1.0;
        m_log_frequency = 5;
        m_method = 0;
        m_gravity << 0, -9.81, 0;
        reset();
    }

    virtual void resetMembers() override {
        setMethod(m_method);

        boundaries.clear();
        grids.clear();
        particles.clear();
        m_objects_ls.clear();

        m_renderVs.clear();
        m_renderFs.clear();
        m_renderNs.clear();

        sandNum = 0;

        if (SCENEMODE != SANDPILE && SCENEMODE != WATER)
            initRigidObjects();
        if (SCENEMODE == CASTLE)
            initParticles_sand_castle();
        else
            initParticles();
        initBoundaries();
        initGrids();

        if (SCENEMODE == CITY)
            initCity();
    }

    virtual void updateRenderGeometry() override {
        for (auto& o : m_objects_ls) {
            o.getMesh(m_renderVs[o.getID()], m_renderFs[o.getID()]);
        }

    }

    virtual bool advance() override;

    virtual void renderRenderGeometry(
        igl::opengl::glfw::Viewer& viewer) override {

        // -------------------------------------
        // Edit here if more rigid body are added
        // -------------------------------------
        int num_meshes = std::max(int(m_objects_ls.size()), 1);
        if (SCENEMODE == CITY)
            num_meshes++;
        viewer.data_list.resize(num_meshes);
        Eigen::MatrixXd color = Eigen::MatrixXd(1, 3);
        color(0, 0) = 1; color(0, 1) = 0; color(0, 2) = 0;
        for (int i = 0; i < m_objects_ls.size(); i++) {
            viewer.data_list[i].show_lines = true;
            viewer.data_list[i].set_face_based(true);
            viewer.data_list[i].point_size = 2.0f;
            viewer.data_list[i].set_mesh(m_renderVs[i], m_renderFs[i]);
            viewer.data_list[i].set_colors(color);
        }
        if (SCENEMODE == CITY) {
            int idx = num_meshes - 1;
            viewer.data_list[idx].show_lines = true;
            viewer.data_list[idx].set_face_based(true);
            viewer.data_list[idx].point_size = 2.0f;
            viewer.data_list[idx].set_mesh(m_renderV_city, m_renderF_city);
        }

        float scale = 3;
        m_particle_matrix = Eigen::MatrixXd(particles.size(), 3);
        m_particle_color_matrix = Eigen::MatrixXd(particles.size(), 3);

        for (int i = 0; i < int(particles.size()); i++)
        {
            m_particle_matrix(i, 0) = particles[i]->x_p.x() / 50;
            m_particle_matrix(i, 1) = particles[i]->x_p.y() / 50;
            m_particle_matrix(i, 2) = particles[i]->x_p.z() / 50;
            m_particle_color_matrix(i, 0) = 1;
            m_particle_color_matrix(i, 1) = 1;
            m_particle_color_matrix(i, 2) = 1;
            if (i < sandNum) {
                m_particle_color_matrix(i, 0) = 0.7078;
                m_particle_color_matrix(i, 1) = 0.6984;
                m_particle_color_matrix(i, 2) = 0.5020;

            }
            else {
                m_particle_color_matrix(i, 0) = 0;
                m_particle_color_matrix(i, 1) = 0.4118;
                m_particle_color_matrix(i, 2) = 0.5804;
            }

        }

        viewer.data_list[0].point_size = scale;
        viewer.data_list[0].points = Eigen::MatrixXd(0, 6);
        viewer.data_list[0].add_points(m_particle_matrix, m_particle_color_matrix);
    }


#pragma region SettersAndGetters
    void updateCannonBall() {
        Eigen::Vector3d momentum;
        momentum << std::cos(m_angle), std::sin(m_angle), 0;
        momentum *= m_force;
        m_objects_ls[0].setLinearVelocity(momentum / m_objects_ls[0].getMass());
    }

    void setAngle(double a) {
        m_angle = a;
    }

    void setForce(double f) {
        m_force = f;
    }

    void setMass(double m) { m_mass = m; }

    void setMethod(int m) {
        m_method = m;
        switch (m_method) {
        case 0:
            m_color = Eigen::RowVector3d(1.0, 0.0, 0.0);
            break;
        case 1:
            m_color = Eigen::RowVector3d(0.0, 1.0, 0.0);
            break;
        case 2:
            m_color = Eigen::RowVector3d(0.0, 0.0, 1.0);
            break;
        default:
            std::cerr << m_method << " is not a valid integrator method."
                << std::endl;
        }

    }

    void setLogFrequency(int f) { m_log_frequency = f; }

    // Initializations
#if Material3D == Sand3D
    void AddParticles() {
        Vector3f v;
        Matrix3f a = Matrix3f::Zero();

        float W_COL, H_COL, P_COL;
        float X_COL, Y_COL, Z_COL;

        int NP = 100;

        if (SCENEMODE == SANDPILE) {
            W_COL = X_GRID / 8.0f;
            H_COL = (Y_GRID - 2 * CUB) * 0.9f;
            P_COL = Z_GRID / 8.0f;

            X_COL = (X_GRID - W_COL) / 2.0f;
            Y_COL = CUB;
            Z_COL = (Z_GRID - P_COL) / 2.0f;
            v.setZero();
            NP = 10;
        }
        else {
            W_COL = X_GRID / 10.0f;
            H_COL = (Y_GRID - 2 * CUB) * 0.3f;
            P_COL = Z_GRID / 10.0f;

            X_COL = (X_GRID - W_COL) / 2.0f;
            Y_COL = CUB + Y_GRID * 0.3;
            Z_COL = Z_GRID - 2 * P_COL - 2 * CUB;
            v = Vector3f(0, -30, 0);
            NP = 100;
        }


        float VOL = W_COL * H_COL * P_COL / static_cast<float>(10000);
        float MASS = VOL * RHO_dry_sand / 100.0f;

        for (int i = 0; i < NP; i++)
        {
            float r1 = ((float)rand() / (RAND_MAX));
            float r2 = ((float)rand() / (RAND_MAX));
            Vector3f pos = Vector3f(r1 * W_COL + X_COL, 1 * H_COL + Y_COL, r2 * P_COL + Z_COL); 
            //particles.push_back(Sand3D(VOL, MASS, pos, v, a));
            particles.push_back(new Sand3D(VOL, MASS, pos, v, a));
        }
        sandNum += NP;

    }

    void AddParticles_water() {
        Vector3f v;
        Matrix3f a = Matrix3f::Zero();

        float W_COL, H_COL, P_COL;
        float X_COL, Y_COL, Z_COL;

        int NP = 100;


        W_COL = X_GRID / 8.0f;
        H_COL = (Y_GRID - 2 * CUB) * 0.9f;
        P_COL = Z_GRID / 8.0f;

        v = Vector3f(-50.0f, 0, 0);
        NP = 10;

        float VOL = W_COL * H_COL * P_COL / static_cast<float>(10000);
        float MASS = VOL * RHO_water / 100.0f;


        H_COL = Y_GRID / 30.f;
        P_COL = Z_GRID / 30.f;

        X_COL = (X_GRID - CUB) * 1.0f;
        Y_COL = (Y_GRID - H_COL) / 2.0f;
        Z_COL = (Z_GRID - P_COL) / 2.0f;

        for (int i = 0; i < NP; i++)
        {
            float r1 = ((float)rand() / (RAND_MAX));
            float r2 = ((float)rand() / (RAND_MAX));
            //Vector3f pos = Vector3f(r1 * W_COL + X_COL, 1 * H_COL + Y_COL, r2 * P_COL + Z_COL);
            Vector3f pos = Vector3f(X_COL, r1 * H_COL + Y_COL, r2 * P_COL + Z_COL);
            //particles.push_back(Sand3D(VOL, MASS, pos, v, a));
            particles.push_back(new Water3D(VOL, MASS, pos, v, a));
        }

    }

    RigidObjectLS* addBall(const double& scale, const Eigen::Vector3d& position, const double& m_mass, const Eigen::Vector3d& velocity, const Eigen::Vector3d& angular_velocity) {
        const std::string path = "sphere.off";
        m_objects_ls.push_back(RigidObjectLS(path));
        auto p_ball = &m_objects_ls.back();

        double radius = scale * scale * 19.7371;
        auto inertia = 0.4 * m_mass * radius * radius * Eigen::Matrix3d::Identity();
        p_ball->setInertia(inertia);
        p_ball->reset();
        p_ball->setScale(scale);
        p_ball->setPosition(position);
        p_ball->setMass(m_mass);
        p_ball->setLinearVelocity(velocity);
        p_ball->setAngularVelocity(angular_velocity);
        p_ball->setRigidObjectType(BALL);
        p_ball->setID(m_objects_ls.size() - 1);

        m_renderVs.push_back(Eigen::MatrixXd());
        m_renderFs.push_back(Eigen::MatrixXi());
        m_renderNs.push_back(Eigen::MatrixXd());
        return p_ball;
    }

    RigidObjectLS* addCube(const double& scale, const Eigen::Vector3d& position, const double& m_mass, const Eigen::Vector3d& velocity, const Eigen::Vector3d& angular_velocity) {
        const std::string path = "cube.off";
        m_objects_ls.push_back(RigidObjectLS(path));
        auto p_cube = &m_objects_ls.back();

        double radius = scale * scale * 1.0;
        p_cube->setInertia((2.0 / 3.0) * m_mass * radius * radius * Eigen::Matrix3d::Identity());
        p_cube->reset();
        p_cube->setScale(scale);
        p_cube->setPosition(position);
        p_cube->setMass(m_mass);
        p_cube->setLinearVelocity(velocity);
        p_cube->setAngularVelocity(angular_velocity);
        p_cube->setRigidObjectType(CUBE);
        p_cube->setID(m_objects_ls.size() - 1);

        m_renderVs.push_back(Eigen::MatrixXd());
        m_renderFs.push_back(Eigen::MatrixXi());
        m_renderNs.push_back(Eigen::MatrixXd());
        return p_cube;
    }

    void initRigidObjects() {
        if (SCENEMODE == CITY) {
            double scale = 0.25;
            auto cube1 = addCube(scale, Eigen::Vector3d(X_GRID / double(100), scale * scale + CUB / double(50), (Z_GRID / 2 + 30) / double(50)), 100.0f, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0));
            auto cube2 = addCube(scale, Eigen::Vector3d(X_GRID / double(100), scale * scale + CUB / double(50), (Z_GRID / 2 + 15) / double(50)), 100.0f, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0));
        }
        else if (SCENEMODE == CASTLE) {
            auto cube1 = addCube(0.2f, Eigen::Vector3d(0, 0, 0), 100000.0f, Eigen::Vector3d(5, 2, 5), Eigen::Vector3d(5, 5, 5));
            // auto cube1 = addCube(0.4f, Eigen::Vector3d(0, 0, 0), 100000.0f, Eigen::Vector3d(5, 2, 5), Eigen::Vector3d(50, 50, 50));
            // auto cube1 = addCube(0.4f, Eigen::Vector3d(0, 0, 0), 100000.0f, Eigen::Vector3d(5, 2, 5), Eigen::Vector3d(5000, 5000, 5000))
        }
        else { // TESTSCENE
            auto ball1 = addBall(0.1f, Eigen::Vector3d(0, 0, 0), 100000.0f, Eigen::Vector3d(5, 5, 5), Eigen::Vector3d(0, 0, 0));
            // auto ball1 = addBall(0.1f, Eigen::Vector3d(0, 0, 0), 1000.0f, Eigen::Vector3d(5, 5, 5), Eigen::Vector3d(0, 0, 0));
            // auto cube1 = addCube(0.4f, Eigen::Vector3d(0, 0, 0), 100000.0f, Eigen::Vector3d(5, 5, 5), Eigen::Vector3d(0, 0, 0));
            // auto cube1 = addCube(0.4f, Eigen::Vector3d(0, 0, 0), 100000.0f, Eigen::Vector3d(5, 5, 5), Eigen::Vector3d(5, 5, 5));
        }
    }

    void initParticles() {
        // -------------------------------
        // Possion sampling
        // -------------------------------
        DefaultPRNG PRNG;
        const auto P_c = GeneratePoissonPoints(50000, PRNG);
        int NP = static_cast <int>(P_c.size());
        float W_COL, H_COL, P_COL;
        float X_COL, Y_COL, Z_COL;
        Vector3f v;

        if (SCENEMODE == CITY) {
            W_COL = X_GRID / 10.0f;
            H_COL = (Y_GRID - 2 * CUB) * 0.3f;
            P_COL = Z_GRID / 10.0f;
            X_COL = (X_GRID - W_COL) / 2.0f;
            Y_COL = CUB + Y_GRID * 0.3;
            Z_COL = Z_GRID - 2 * P_COL - 2 * CUB;
            v = Vector3f(0.0f, -30.0f, 0.0f);  						// Initial velocity
        }
        else { // TESTSCENE
            W_COL = X_GRID / 8.0f;
            H_COL = (Y_GRID - 2 * CUB) * 0.9f;
            P_COL = Z_GRID / 8.0f;
            X_COL = (X_GRID - W_COL) / 2.0f;
            Y_COL = CUB;
            Z_COL = (Z_GRID - P_COL) / 2.0f;
            v = Vector3f::Zero();
        }

        float VOL = W_COL * H_COL * P_COL / static_cast<float>(10000);
        float MASS = VOL * RHO_dry_sand / 100.0f;
        Matrix3f a = Matrix3f::Zero();

        for (int p = 0; p < NP; p++)
        {
            float r = ((float)rand() / (RAND_MAX));
            Vector3f pos = Vector3f(P_c[p].x * W_COL + X_COL, P_c[p].y * H_COL + Y_COL, r * P_COL + Z_COL); 
            //particles.push_back(Sand3D(VOL, MASS, pos, v, a));
            if (SCENEMODE != WATER)
                particles.push_back(new Sand3D(VOL, MASS, pos, v, a));
            else
                particles.push_back(new Sand3D(VOL, MASS * RHO_water /RHO_dry_sand, pos, v, a));
        }
        if (SCENEMODE != WATER)
            sandNum += NP;
    }

    void initParticles_sand_castle() {
        // -------------------------------
        // Initialize sand castle
        // -------------------------------
        std::string path = "model/castle.stl";

        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::MatrixXd N;
        if (!igl::readSTL(path, V, F, N)) {
            path = "../" + path;
            if (!igl::readSTL(path, V, F, N)) {
                path = "../" + path;
                igl::readSTL(path, V, F, N);
            }
        }

        vector<double> maxCoor(V.cols(), numeric_limits<int>::min());
        vector<double> minCoor(V.cols(), numeric_limits<int>::max());

        int range = 100;
        for (int i = 0; i < V.rows(); i++) {
            for (int j = 0; j < V.cols(); j++) {
                V(i, j) = range * V(i, j);
                if (maxCoor[j] < V(i, j)) {
                    maxCoor[j] = V(i, j);
                }
                if (minCoor[j] > V(i, j)) {
                    minCoor[j] = V(i, j);
                }
            }
        }

        /*
        std::cout << "V.rows(): " << V.rows() << ", F.rows(): " << F.rows() << std::endl;
        std::cout << "V.cols(): " << V.cols() << ", F.cols(): " << F.cols() << std::endl;
        std::cout << "maxCoor: ";
        for (int i = 0; i < V.cols(); i++) {
            std::cout << maxCoor[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "minCoor: ";
        for (int i = 0; i < V.cols(); i++) {
            std::cout << minCoor[i] << " ";
        }
        std::cout << std::endl;
        */

        int scale = 100;
        vector<double> inner(int(maxCoor[1] - minCoor[1]) / scale + 1, 0);
        vector<vector<double>> maxZ(int(maxCoor[0] - minCoor[0]) / scale + 1, inner);
        for (int i = 0; i < V.rows(); i++) {
            int x = int(V(i, 0) - minCoor[0]) / scale;
            int y = int(V(i, 1) - minCoor[1]) / scale;
            double z = (V(i, 2) - minCoor[2]) / scale;
            maxZ[x][y] = max(z, maxZ[x][y]);
        }

        int NP = 0;
        float VOL = 1;
        float MASS = VOL * RHO_dry_sand / 100.0f;

        Vector3f v = Vector3f::Zero();							// Initial velocity
        Matrix3f a = Matrix3f::Zero();

        float vol_size = 1;
        int perLayerNum = 3000;

        DefaultPRNG PRNG;
        int X_offset = (X_GRID - (maxCoor[0] - minCoor[0]) / scale) / 2.0f;
        int Z_offset = (Z_GRID - (maxCoor[1] - minCoor[1]) / scale) / 2.0f;
        for (int z = 1; z < (maxCoor[2] - minCoor[2]) / scale + 1; z++) {
            const auto P_c = GeneratePoissonPoints(perLayerNum, PRNG);

            int num = static_cast <int>(P_c.size());
            for (int i = 0; i < num; i++) {
                int x = (P_c[i].x * (maxCoor[0] - minCoor[0])) / scale;
                int y = (P_c[i].y * (maxCoor[1] - minCoor[1])) / scale;

                if (z <= maxZ[x][y]) {
                    float r = ((float)rand() / (RAND_MAX));
                    Vector3f pos = Vector3f((P_c[i].x * (maxCoor[0] - minCoor[0])) / scale + X_offset,
                        max(double(CUB), z * vol_size + CUB + r - 0.5),
                        (P_c[i].y * (maxCoor[1] - minCoor[1])) / scale + Z_offset); // notice that the orientation of the object file is different from the setting in this project
                    //particles.push_back(Sand3D(VOL, MASS, pos, v, a));
                    particles.push_back(new Sand3D(VOL, MASS, pos, v, a));
                    NP++;
                }
            }
        }
        sandNum += NP;
    }

    void initCity() {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::MatrixXd N_faces;

        std::string path = "model/city.obj";
        if (!igl::readOBJ(path, V, F)) {
            path = "../" + path;
            if (!igl::readOBJ(path, V, F)) {
                path = "../" + path;
                igl::readOBJ(path, V, F);
            }
        }

        double scale = 7000;
        vector<double> maxCoor(V.cols(), numeric_limits<double>::min());
        vector<double> minCoor(V.cols(), numeric_limits<double>::max());
        for (int i = 0; i < V.rows(); i++) {
            for (int j = 0; j < V.cols(); j++) {
                V(i, j) = V(i, j) * scale;
                if (maxCoor[j] < V(i, j)) {
                    maxCoor[j] = V(i, j);
                }
                if (minCoor[j] > V(i, j)) {
                    minCoor[j] = V(i, j);
                }
            }
        }
        float xScale = (X_GRID - CUB * 2) / (maxCoor[0] - minCoor[0]);
        float zScale = (Z_GRID - CUB * 2) / (maxCoor[2] - minCoor[2]);

        /*
        std::cout << "V.rows(): " << V.rows() << ", F.rows(): " << F.rows() << std::endl;
        std::cout << "V.cols(): " << V.cols() << ", F.cols() " << F.cols() << std::endl;

        std::cout << "maxCoor: ";
        for (int i = 0; i < V.cols(); i++) {
            std::cout << maxCoor[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "minCoor: ";
        for (int i = 0; i < V.cols(); i++) {
            std::cout << minCoor[i] << " ";
        }
        std::cout << std::endl;
        */

        for (int i = 0; i < V.rows(); i++) {
            V(i, 0) = (-V(i, 0) + maxCoor[0]) * xScale + CUB;
            V(i, 1) = V(i, 1) - minCoor[1] + CUB;
            V(i, 2) = (-V(i, 2) + maxCoor[2]) * zScale + CUB;
            
            // ignore the little bump up
            if (V(i, 1) < 5) {
                V(i, 1) = CUB;
            }
        }

        igl::per_face_normals(V, F, N_faces);

        // -------------------------------
        // version 1: consider city as boundary
        // -------------------------------
        //for (int i = 0; i < F.rows(); i++) {
        //    std::vector<Vector3f> Corners;
        //    //for (int j = 0; j < F.cols(); j++) {
        //    if (1 - N_faces(i, 1) > eps) continue;
        //    for (int j = 0; j < 3; j++) {
        //        int idx = F(i, j);
        //        Corners.push_back(Vector3f(V(idx, 0), V(idx, 1), V(idx, 2)));
        //    }

        //    Vector3f faceNormal = Vector3f(N_faces(i, 0), N_faces(i, 1), N_faces(i, 2));
        //    boundaries.push_back(Boundary3D((faceNormal), 3, Corners));
        //    bool status = true;

        //    for (int q = 0; q < 3; q++) {
        //        status = status || abs(abs(faceNormal[q])-1) < eps;
        //    }
        //    if (!status) {
        //        std::cout << "!status: " << N_faces(i, 0) << " " << N_faces(i, 1) << " " << N_faces(i, 2) << std::endl;
        //        std::cout << std::endl;
        //    }
        //}

        // -------------------------------
        // version 2: develop the level set for the city
        // -------------------------------
        vector<pair<int, tuple<int, int, int, int>>> upperFaces;
        vector<vector<int>> maxY(X_GRID + 1, vector<int>(Z_GRID + 1, 0));
        for (int i = 0; i < F.rows(); i++) {
            if (1 - N_faces(i, 1) > eps) continue;
            int x_min = numeric_limits<int>::max();
            int x_max = numeric_limits<int>::min();
            int z_min = numeric_limits<int>::max();
            int z_max = numeric_limits<int>::min();
            for (int j = 0; j < 3; j++) {
                int idx = F(i, j);
                x_min = min(x_min, int(V(idx, 0)));
                x_max = max(x_max, int(V(idx, 0)));
                z_min = min(z_min, int(V(idx, 2)));
                z_max = max(z_max, int(V(idx, 2)));
            }
            for (int x = x_min; x <= x_max; x++) {
                for (int z = z_min; z <= z_max; z++) {
                    maxY[x][z] = max(maxY[x][z], int(V(F(i, 0), 1)));
                }
            }
            if (x_max == x_min || z_max == z_min) continue;
            upperFaces.push_back(make_pair(i, make_tuple(x_min, x_max, z_min, z_max)));
        }

        normalField = vector<Vector3f>((X_GRID + 1) * (Y_GRID + 1) * (Z_GRID + 1), Vector3f(0, 0, 0));
        valueField = vector<float>((X_GRID + 1) * (Y_GRID + 1) * (Z_GRID + 1), 1);
        for (size_t i = 0; i < upperFaces.size(); i++) {
            int idx_f = upperFaces[i].first;
            int y_f = int(V(F(idx_f, 0), 1));
            if (y_f <= CUB) continue;

            int x_min = get<0>(upperFaces[i].second);
            int x_max = get<1>(upperFaces[i].second);
            int z_min = get<2>(upperFaces[i].second);
            int z_max = get<3>(upperFaces[i].second);
            for (int x = x_min; x <= x_max; x++) {
                for (int z = z_min; z <= z_max; z++) {
                    if (y_f < maxY[x][z]) continue; // only consider the upper most face
                    for (int y = 0; y <= y_f; y++) {
                        int idx = (X_GRID + 1) * (Y_GRID + 1) * z + (X_GRID + 1) * y + x;
                        valueField[idx] = -1;
                    }
                }
            }
        }

        for (int x = 0; x < X_GRID + 1; x++) {
            for (int y = 0; y < Y_GRID + 1; y++) {
                for (int z = 0; z < Z_GRID + 1; z++) {
                    int idx = (X_GRID + 1) * (Y_GRID + 1) * z + (X_GRID + 1) * y + x;
                    if (valueField[idx] > 0) continue;
                    int x_count_posi = 1;
                    int x_count_neg = 1;
                    int z_count_posi = 1;
                    int z_count_neg = 1;
                    int y_count_posi = 1;

                    while (valueField[idx + x_count_posi] <= 0) { x_count_posi++; }
                    while (valueField[idx - x_count_neg] <= 0) { x_count_neg++; }
                    while (valueField[idx + z_count_posi * (X_GRID + 1) * (Y_GRID + 1)] <= 0) { z_count_posi++; }
                    while (valueField[idx - z_count_neg * (X_GRID + 1) * (Y_GRID + 1)] <= 0) { z_count_neg++; }
                    while (valueField[idx + y_count_posi * (X_GRID + 1)] <= 0) { y_count_posi++; }

                    if (x_count_posi <= min(min(min(x_count_neg, z_count_posi), z_count_neg), y_count_posi)) {
                        normalField[idx] = Vector3f(1, 0, 0);
                        valueField[idx] = -x_count_posi;
                    }
                    else if (x_count_neg <= min(min(min(x_count_posi, z_count_posi), z_count_neg), y_count_posi)) {
                        normalField[idx] = Vector3f(-1, 0, 0);
                        valueField[idx] = -x_count_neg;
                    }
                    else if (z_count_posi <= min(min(min(x_count_posi, x_count_neg), z_count_neg), y_count_posi)) {
                        normalField[idx] = Vector3f(0, 0, 1);
                        valueField[idx] = -z_count_posi;
                    }
                    else if (z_count_neg <= min(min(min(x_count_posi, x_count_neg), z_count_posi), y_count_posi)) {
                        normalField[idx] = Vector3f(0, 0, -1);
                        valueField[idx] = -z_count_neg;
                    }
                    else {
                        normalField[idx] = Vector3f(0, 1, 0);
                        valueField[idx] = -y_count_posi;
                    }
                }
            }
        }

        for (int x = 0; x < X_GRID + 1; x++) {
            for (int y = 0; y < CUB; y++) {
                for (int z = 0; z < Z_GRID + 1; z++) {
                    int idx = (X_GRID + 1) * (Y_GRID + 1) * z + (X_GRID + 1) * y + x;
                    if (valueField[idx] > 0) {
                        valueField[idx] = -1;
                        normalField[idx] = Vector3f(0, 1, 0);
                    }
                }
            }
        }

        /*
        // Visualization of normal/value field
        int y = 20;
        for (int z = 0; z < Z_GRID + 1; z++) {
            for (int x = 0; x < X_GRID + 1; x++) {
                int idx = (X_GRID + 1) * (Y_GRID + 1) * z + (X_GRID + 1) * y + x;
                if (valueField[idx] > 0)
                    std::cout << "-";
                else {
                    if(normalField[idx][0] > 1-eps)
                        std::cout << ">";
                    else if (-normalField[idx][0] > 1 - eps)
                        std::cout << "<";
                    else if (normalField[idx][2] > 1 - eps)
                        std::cout << "+";
                    else if (-normalField[idx][2] > 1 - eps)
                        std::cout << "^";
                    else
                        std::cout << ".";
                }
            }
            std::cout << std::endl;
        }
        */

        // adapt the scale accordingly
        V = V / 50;

        m_renderV_city = V;
        m_renderF_city = F;
    }

#endif

    void initBoundaries() {
        std::vector<Vector3f> Corners;
        /* x = 0 */
        Corners.push_back(Vector3f(CUB, CUB, CUB));
        Corners.push_back(Vector3f(CUB, Y_GRID - CUB, CUB));
        Corners.push_back(Vector3f(CUB, Y_GRID - CUB, Z_GRID - CUB));
        Corners.push_back(Vector3f(CUB, CUB, Z_GRID - CUB));
        boundaries.push_back(Boundary3D(Vector3f(1, 0, 0), 2, Corners));
        Corners.clear();

        /* y = 0 */
        Corners.push_back(Vector3f(CUB, CUB, CUB));
        Corners.push_back(Vector3f(CUB, CUB, Z_GRID - CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, CUB, Z_GRID - CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, CUB, CUB));
        boundaries.push_back(Boundary3D(Vector3f(0, 1, 0), 2, Corners));
        Corners.clear();

        /* z = 0 */
        Corners.push_back(Vector3f(CUB, CUB, CUB));
        Corners.push_back(Vector3f(CUB, Y_GRID - CUB, CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, Y_GRID - CUB, CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, CUB, CUB));
        boundaries.push_back(Boundary3D(Vector3f(0, 0, 1), 2, Corners));
        Corners.clear();

        /* x =  X_GRID*/
        Corners.push_back(Vector3f(X_GRID - CUB, CUB, CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, Y_GRID - CUB, CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, Y_GRID - CUB, Z_GRID - CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, CUB, Z_GRID - CUB));
        boundaries.push_back(Boundary3D(Vector3f(-1, 0, 0), 2, Corners));
        Corners.clear();

        /* y = Y_GRID */
        Corners.push_back(Vector3f(CUB, Y_GRID - CUB, CUB));
        Corners.push_back(Vector3f(CUB, Y_GRID - CUB, Z_GRID - CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, Y_GRID - CUB, Z_GRID - CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, Y_GRID - CUB, CUB));
        boundaries.push_back(Boundary3D(Vector3f(0, -1, 0), 2, Corners));
        Corners.clear();

        /* z = Z_GRID */
        Corners.push_back(Vector3f(CUB, CUB, Z_GRID - CUB));
        Corners.push_back(Vector3f(CUB, Y_GRID - CUB, Z_GRID - CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, Y_GRID - CUB, Z_GRID - CUB));
        Corners.push_back(Vector3f(X_GRID - CUB, CUB, Z_GRID - CUB));
        boundaries.push_back(Boundary3D(Vector3f(0, 0, -1), 2, Corners));
        Corners.clear();
    }


    void initGrids() {
        for (int z = 0; z <= Z_GRID; z++)
            for (int y = 0; y <= Y_GRID; y++)
                for (int x = 0; x <= X_GRID; x++)
                    grids.push_back(Grid3d(Vector3f((float)x, (float)y, (float)z)));
    }


#pragma endregion SettersAndGetters

private:
    int m_method;  // id of integrator to be used (0: analytical, 1: explicit
                   // euler, 2: semi-implicit euler)
    double m_angle;
    double m_force;
    double m_mass;

    Eigen::Vector3d m_gravity;
    std::vector<Eigen::MatrixXd> m_renderVs;
    std::vector<Eigen::MatrixXi> m_renderFs;
    std::vector<Eigen::MatrixXd> m_renderNs;

    Eigen::MatrixXd m_renderV_city;  // vertex positions for rendering city
    Eigen::MatrixXi m_renderF_city;  // face indices for rendering city
    Eigen::MatrixXd m_renderN_city;  // Normal indices for rendering city

    Eigen::MatrixXd m_particle_matrix;
    Eigen::MatrixXd m_particle_color_matrix;

    int m_log_frequency;  // how often should we log the COM in the GUI
    Eigen::RowVector3d m_color;

    std::vector<RigidObjectLS> m_objects_ls;
    std::vector<Boundary3D> boundaries;
    std::vector<Grid3d> grids;
    std::vector<int> valid_grid_indexes; // trivial speed up to get openmp work for grids
    //std::vector<Material3D> particles;
    std::vector<MPMParticle3D*> particles;

    std::vector<Vector3f> normalField;
    std::vector<float> valueField;

    int sandNum = 0;


    // member functions
    bool advance_mpm();
    bool advance_rigid_bodies();
    bool ENABLE_TIMER = false;

    void P2G();										// Transfer from Particles to Grid nodes
    void GetValidGrids();                           // Speed up valid grid acquisition
    void UpdateGrids();
    void G2P();										// Transfer from Grid nodes to Particles
    void UpdateParticles();
    void ResetGrid();

    void cityCollision();
    void cityFriction();
};

