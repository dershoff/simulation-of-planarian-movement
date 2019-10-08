# Project Description
This repository is a supplementary material for the publication on scientific research (please refer to it for more details on the scientific part of the project[, link](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3275291))
### "Emergence of Bilateral Symmetry from Chiral Components in the Planarian Epidermis" 

This part of the research attempts to explain movement patterns of Planaria Worms using bottom-up approach: knowing some details on cilia beating, we predict worm movement trajectory. Below we provide a very short recap on how the simulation is done; for a detailed explanation please read __/docs/Simulation.pdf__.  

The simulation views the worm as a solid body that translates and rotates as a whole without any deformation. The movement is due to many cilia that beat with a user-specified coherency; each cilia beats with a certain force, direction, and on/off state. In each iteration step the total translational force and rotational force (torque) is calculated form a newly drawn configuration of cilia force parameters.
These forces re then translated into movement of the worm and the trajectory is displayed.

A user can set the angular bias of the cilia  beating  for a chosen subpart of the worm, which will most likely result in asymmetry in the worm trajectory. In this reaserarch we use an asymmetric directional bias that has been obtained experimentally (using fluorescent miroscopy imaging) and predict the movement trajectories;  the results are well consistent (qualitatively) with the real trajectories, and can even predict movement patterns that later were confirmed experimentally.

# Installation
No installation required; just open it in your matlab session, or add it to Matlab path (right click on the folder and choose "Add to Path").

# Usage
1. First create settings for the simulation.
    ``` Matlab 
    opts.simulation.frames = 1000;
    opts.simulation.showFrames = 0;
    opts.simulation.saveVideo = 0;

    opts.worm.length = 8; %mm
    opts.worm.bodyRatio = 2.8; %length to width, normal value = 2.8
    opts.worm.initialPosition = [5, 3, 0];
    opts.worm.orientation = 90; %90 = upwards

    opts.worm.strokeForce = 0.2; %arbitrary unit
    opts.worm.mass = 10; %arbitrary
    opts.worm.dt = 1; %time step for integration
    ```
2. Create a worm:
    ``` Matlab 
    w = gen_worm(opts);
    w.singleShow; %visualise

    ```
3. Create settings for simulation:
    ``` Matlab 
    opts.simulation.frames = 100;
    opts.simulation.showFrames = 1;
    opts.simulation.showStamps = 0;
    opts.simulation.saveVideo = 0;
    opts.simulation.VideoName = 'temp_test';
    ```

3. Run simulation. Iterations are done inside the funciton "runStuff"; the acrtual workhorse is the method of the worm class "worm.move(), which updates the coordinates due to cilia forces and updates the cilia force parameters.
    ``` Matlab 
    plotOn =  1; % to visualize iteraitno steps
    runStuff(w, opts, plotOn);
    ````
4. The output parameters can be accumulated and then plotted.
    ``` Matlab 
    TODO;
    ```

Look for more exmaples in the coresponding folder.

Feel free to use the code for your own purposes; please dont forget to reference the publication in the link above.

