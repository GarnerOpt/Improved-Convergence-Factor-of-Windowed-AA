%==========================================================================
% Description: 
% This script runs all of the numerical experiments conducted
% with the TME model. This performs all of the numerical WITHOUT loading any 
% the numerical results to reduce computational time. 
%
% Note, the run times will vary from machine to machine. So, the box-plots
% will not look exactly the same. 
% 
% To generate the figures in the paper run: 
%     Nonlinear_Operation_Simulations_SHORT.m
%==========================================================================
clc; clear all; 

close all;
NonLinear_TME_Operator_Experiment_1

clear all; 
NonLinear_TME_Operator_Experiment_2

clear all; 
NonLinear_TME_Operator_Experiment_3

clear all; 
NonLinear_TME_Operator_Experiment_4

