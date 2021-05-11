# ProductionOptRig
Supporting codes of the paper: Implementation of Steady-state Real-time Optimization Using Transient Measurements on Experimental Rig

The files here are divided in two main folders: 
1. Experimental Data: The files containt the data used for generating the plots in the paper. 
2. MockUpLoop: 

In this file, we have a mockup loop representing the experimental rig, where a dynamic model is used for representing the system. The files 

a) InitializationLabViewMain.m | b) LabViewMain.m 

are the same files, with the same tuning parameters that are used in the actual experimental rig. Therefore, the implementation of the three methods discussed in the paper (SSRTO, DRTO, and ROPA) can be studied by running these files. For obtaining the simulation results, simply run the file Main.m included in the folder with the same name as the method of interest.

N.B.: All the data presented in the paper is experimental, no simulational results were used in the paper analysis.
