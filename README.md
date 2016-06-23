# SMC-EM

This is a project implementing the Sequential Monte Carlo - Expectation Maximization 
(SMC-EM) algorithm described in T. T. Ashley and S. B. Andersson, “Method for simultaneous
localization and parameter estimation in particle tracking experiments,” Physical Review 
E, vol. 92, no. 5, pp. 052707–1–19, Nov. 2015.

DEPENDENCIES:
Requires the boost and tiff modules.

USAGE:
Use requires a bit of familiarity with C++ but not too much. In short, make modifications 
to the beginning part of main.cpp to define your optimization parameters, your optical
parameters, your desired motion and observation models, and set initial conditions. Then
build the project, load the needed modules, and execute the built program.

Execution requires paths to the input data, trajectory data and where you want the output 
to go. Thus, call it as something like
./main $DATA_DIR $RESULTS_DIR $SENSOR_PATH

INPUTS
The data. This can be a sequence of images (wide-field imaging) or time-series (single-
pixel/confocal). The files should be in $DATA_DIR.

The sensor trajectories. This is a data file specifying where the data was taken. In the
wide-field case, it would be the position of the center of each pixel. 

OUTPUTS
Three different files are made at each EM iteration. locxx is a file containing the 
discrete approximation to the PDF of the particle position at each time step. paramxx 
contains the estimated parameters; contents depends on the different modules chosen. 
jointxx contains the joint densities.

Code was originally written by Trevor T. Ashley.