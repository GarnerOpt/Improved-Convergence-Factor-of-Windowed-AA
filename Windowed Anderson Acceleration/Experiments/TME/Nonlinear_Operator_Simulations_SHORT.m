%==========================================================================
% Description: 
%-------------------------------------------------------------------------
% This script runs all of the numerical experiments conducted
% with the TME model and generates all of the associated figures in the
% paper. 
%
% For the experiments which take longer the data associated with the
% solutions are imported and the figures made. If the entire computations
% is desired to be performed run Nonlinear_Operator_Simulations.m
%==========================================================================

clc; clear all; close all;
NonLinear_TME_Operator_Experiment_1

clear all;
NonLinear_TME_Operator_Experiment_2_Only_Plots

clear all;
NonLinear_TME_Operator_Experiment_3_Only_Plots

clear all;
NonLinear_TME_Operator_Experiment_4_Only_Plots

