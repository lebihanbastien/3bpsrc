%--------------------------------------------------------------------------
% Init script for all examples 
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------

%% Reboot
clear all;
close all;


%% Add subfolders to the path
addpath('./computation');
addpath('./richardson');
addpath('./data');
addpath('./init');
addpath('./ode');
addpath('./plot');
addpath('./results');
addpath('./scripts');
addpath('./diffcorr');

%% Data loading (abacuses)
load halo_init_matrix_EML2 halo_init_EML2;
load halo_init_matrix_EML1 halo_init_EML1;

%% Constants init
cst = constants_init();

%% Parameters init (to default values, see values within routine)
default = parameters_default_init(cst);