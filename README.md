# Ammonia-RO-model


My original code is pretty messy, chatGPT helped to clean the code and reorganize the notation. But I haven't test the cleaned version. It should work...


The file 'NH4Cl_varing_pH.mlx' is the simplest version considering only pure NH4Cl solution with complete loops for concentration polarization.
This file is a good start to play with. Cl- partitioning coefficient is simplified as a fixed value here


File 'NH4Cl_varying_TAN_Concentration.mlx' is the code for calculating Fig S2 and S3. It involves the varying Cl- partitioning coefficient as a funtion of total
cationic charge concentration.


File 'Syn_UF_effluent.mlx' considered more complex feed matrices containing more competing ions.


Folder 'Fitting' contains the script for searching for optimal fittted parameters. Run 'fit_separate.m' it will automatically call other functions in this folder. 
This code only optimize parameters based on the pure NH4Cl varying pH condition, so the output provide a preliminary clue for the parameters, which were further manully adjusted to fit
the varing pressure and varying TAN concentration conditions.

Folder 'Sobol SA' contains the script of doing sensitivity analysis. Run 'sa_sobol.m'
