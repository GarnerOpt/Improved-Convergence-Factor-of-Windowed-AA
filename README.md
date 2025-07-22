# Improved-Convergence-Factor-of-Windowed-AA
This repository stores the code for the numerical experiments in our paper - "Improved Convergence Factor of Windowed Anderson Acceleration for Symmetric Fixed-Point Iterations." 


-------------
Description: 
-------------
This README file explains how to reproduce the results contained in the paper, "Improved Convergence Rates of Anderson Acceleration for a Large
Class of Fixed-Point Iterations". All code is written in MATLAB R2021b and running the code requires access to a MATLAB license. 

-------
Steps:
-------
(1) Download the code in the folder called, Windowed Anderson Acceleration. 

(2) Before running the necessary scripts make sure the folder and all its subfolders are on the the file path. 

(3) Run the following scripts from the command line to reproduce the figures in the paper, or similar figures under the same data settings: 
	* Linear_Operator_Simulations.m
	* Nonlinear_Operator_Simulations.m
	* Nonlinear_Operator_Simulations_SHORT.m  (recommended, see (2) and (3) below)

-----------------	
IMPORTANT NOTES:
-----------------	
(1) Each of these .m files produces a different set of numerical results. The first, Linear_Operator_Simulations.m, produces all of the experiments with linear operators. 
More precisly, running this script reproduces Figures 2 and 3 in the paper. 

(2) Nonlinear_Operator_Simulations.m and Nonlinear_Operator_Simulations_SHORT.m  both reproduce all of the experiments on the TME model. 
The first, Nonlinear_Operator_Simulations.m, runs all of the experiments in their entirety. For this reason, it takes about an hour to run.

(3) Nonlinear_Operator_Simulattions_SHORT.m presents the figures presented in the paper by loading the data computed from these experiments. 
Thus, this script takes significantly less time to run when compared to Nonlinear_Operator_Simulations.m
