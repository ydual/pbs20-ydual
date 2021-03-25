# PBS 2020 - Course Project Group 15 - Sand Simulation
Developed by Yinwei Du, Shaohui Liu, and Linfei Pan.

# Project Overview
In this repository, we provide our implementation for our MPM solver based sand simulation. Our codebase is developed directly upon the PBS excercise framework, so the installation step is basically the same. In the demo, the time step is set to be relatively small (only 10ms) to achieve robust and high quality results. However, it is possible to enlarge the time step so that our implementation can run in real time.

# Run demo
After building the project, please run `./build/project` at the project folder. We provide three demo scenes here: 1) Sand flooding into a city. 2) Sand crushing into a castle. 3) Basic test scene (cannon ball crushing into a sand tower). 4) Sand constant droping and piling up. You can run the corresponding demo with its ID. 

```
./build/project ${DEMO_ID}
```

# Final Presentation
Some demo videos can be found via [here](https://ethz-my.sharepoint.com/:p:/g/personal/duyin_ethz_ch/ETgkLMrzCR9MqU-RM0McmuUBOeZEJv-82slQR29n9fm2Cg?e=jghtAr).
Note: You may need ETH account to view the shared file.

# Milestones
* [x] Material Point Method framework in 2D/3D.
* [x] Drucker-Prager elastoplasticity model.
* [x] Coupled with Rigid Bodies (collision detection, momentum and angular momentum updates)
* [x] Parallelization with OpenMP.
* [x] Sand castle and large city Settled as test scenes.
* [x] Collision detection with static scenes (sand / rigid bodies)
* [x] Triangular meshes of the city transferred to signed distance fields to speed up collision handling.
* [x] Setting up the whole scene and parameter tuning.
* [x] Basic rendering with blender.

# Known bugs for visualization
* We found that the GUI from the PBS excercise framework is sometimes unstable so we kindly ask you to carefully operate your mouse over the GUI. That is to say, do not move your mouse so quickly, and do not minimize the window during simulation. 
* Sand particles sometimes can get stick on the boundaries on the side and top


# Acknowledgement
We thank all the lecturers and TAs for providing helpful advice on choosing MPM over SPH as our base method and tips for debugging.

