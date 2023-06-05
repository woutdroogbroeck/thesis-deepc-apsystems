# thesis-deepc-apsystems
Control-oriented Type 1 Diabetes & Artificial Pancreas Simulator for matlab

Thesis topic: Data-Enabled Predictive Control for Artificial Pancreas Systems

## Code structure

The folder [controllers](controllers/) contains the classes of the implemented controllers: LQR, MPC and DeePC (with and without integral action). In [controllers/src](controllers/src/) are some functions that the controllers use (e.g. creating hankel matrix from offline input/output data), as well as a class for state estimation with the Kalman Filter. 

[patient](patient/) contains the file for the Patient class, which can be used to define a virtual patient object with either the minimal model or uva/padova model. [patient/parameters](patient/parameters/) contains the .mat parameter files for the minimal and maximal patients. 

[pump](pump/) and [sensor](sensor/) contain the classes to add a CGM and insulin pump to the patient, as well as their parameter files. 


## Run the code

[setup_path](setup_path.m) adds all the classes, functions from the folders and subfolders to the matlab path.

The controllers can be initialized in [main_controllers.m](main_controllers.m), defining a control policy that takes as input the current state/output and meal announcement object. 

In [main.m](main.m) the simulation parameters can be set such as the integration step, which model to use, linear or nonlinear simulation... The patient_id specifies which patient from the parameter files in [patient/parameters](patient/parameters/) is simulated. For the minimal model there are 3 patients. For uvapadova model there are 30 patients, id 1-10 adolescents, 11-20 adults and 21-30 children. id 31, 32 and 33 are patients with the average parameters of the adolescents, adults and children respectively. This patient can be loaded in using [load_patient](patient/src/load_patient.m) 

Sensor and pump can be added for realistic simulations by commenting out the lines to add them to the patient. 

meals is an array containing the meal sizes for the maximal model. and the fisher disturbance for the minimal model. 

meal_announcement is an array of structs that can be specified as zero if now meal announcement is wanted. Or given a struct with meal.size and meal.minutes containing information about the meal size and minutes until the meal will be consumed (as the commented out example). If the meal is announced 30 minutes before, it should be at index idx_meal-30/dt in meal_announcement. 

simulate in [main.m](main.m) is a function handle of [src/simulation.m](src/simulation.m) that performs a simulation with a given patient and control policy over the number of simulation steps specified by the begin and end time and integration step. This function returns a result struct that contains the usefull variables for plotting
