%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, 2010-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Basic Multiobjective Optimization Example Case
%
%           Aero Goal 1: max(expected(L)/expected(D)) with free transition
%           Aero Goal 2: max(expected(L)/expected(D)) with trip transition
%
%           Structural constraint: t/c = 0.24 (Inwind tip)
%
%           Parallel Execution Enabled
%
%           Drives other cases
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all; clc;

% For sigma 20 (robust design)
N0_new_prob_case_finalrun_sigma20

% For sigma 00 (point design)
N0_new_prob_case_finalrun_sigma00


    



