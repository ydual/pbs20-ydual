#include <igl/writeOFF.h>
#include "MPMSim3D.h" 
#include "Gui.h"

/*
 * GUI for the cannonball simulation. This time we need additional paramters,
 * e.g. which integrator to use for the simulation and the force applied to the
 * cannonball, and we also add some more visualizations (trajectories).
 */
class MPMGui : public Gui {
   public:
    // simulation parameters
    float m_angle = 1.047f;
    float m_force = 30.0f;
    float m_dt = DT;
    float m_mass = 1.0;
    int m_log_frequency = 30;

    MPMSim3D *p_MPMSim = NULL;

    const vector<char const *> m_integrators = {"Analytic", "Explicit Euler", "Symplectic Euler"};
    int m_selected_integrator = 0;

    MPMGui(const int &scene_mode) {
        // create a new cannonball simulation, set it in the GUI,
        // and start the GUI
        p_MPMSim = new MPMSim3D(scene_mode);
        setSimulation(p_MPMSim);

        // show vertex velocity instead of normal
        callback_clicked_vertex = [&](int clickedVertexIndex,
                                      int clickedObjectIndex,
                                      Eigen::Vector3d &pos,
                                      Eigen::Vector3d &dir) {
            RigidObject &o = p_MPMSim->getObjects()[clickedObjectIndex];
            pos = o.getVertexPosition(clickedVertexIndex);
            dir = o.getVelocity(pos);
        };
        start();
    }

    virtual void updateSimulationParameters() override {
        // change all parameters of the simulation to the values that are set in
        // the GUI
        p_MPMSim->setForce(m_force);
        p_MPMSim->setAngle(m_angle);
        p_MPMSim->setTimestep(m_dt);
        p_MPMSim->setMass(m_mass);
        p_MPMSim->setMethod(m_selected_integrator);
        p_MPMSim->setLogFrequency(m_log_frequency);
    }

    virtual void clearSimulation() override {}

    virtual bool childKeyCallback(igl::opengl::glfw::Viewer &viewer,
                                  unsigned int key, int modifiers) override {
        switch (key) {
            case 'e':
            case 'E':
                return true;
            // cicle through different integrators
            case '>':
                m_selected_integrator++;
                m_selected_integrator %= m_integrators.size();
                return true;
            case '<':
                m_selected_integrator--;
                m_selected_integrator =
                    (m_integrators.size() + m_selected_integrator) %
                    m_integrators.size();
                return true;
        }
        return false;
    }

    virtual void drawSimulationParameterMenu() override {
        ImGui::SliderAngle("Angle", &m_angle, -180.0f, 180.0f);
        ImGui::InputFloat("Force", &m_force, 0, 0);
        ImGui::InputFloat("Mass", &m_mass, 0, 0);
        ImGui::InputFloat("dt", &m_dt, 0, 0);
        ImGui::Combo("Integrator", &m_selected_integrator, m_integrators.data(),
                     m_integrators.size());
        ImGui::InputInt("Log Frequency", &m_log_frequency, 0, 0);
    }
};

int main(int argc, char *argv[]) {
    // get scene mode
    int scene_mode = 0;
    if (argc >= 2) {
        std::string k = argv[1];
        std::istringstream ss(k);
        ss >> scene_mode;
    }

    // create a new instance of the GUI for the MPM simulation
    new MPMGui(scene_mode);
    return 0;
}

