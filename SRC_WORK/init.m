%--------------------------------------------------------------------------
% Init script for all computation. Should be called at the beginning of
% each computation.
%
% BLB 2016
%--------------------------------------------------------------------------

%% Reboot
clear all;
close all;


%% Add subfolders to the path
% Recursively add all subfolders
addpath(genpath('.'));

%% Data loading (abacuses)
load halo_init_matrix_EML2 halo_init_EML2;
load halo_init_matrix_EML1 halo_init_EML1;
load nro_init_EML1 nro_init_EML1;
load nro_init_EML2 nro_init_EML2;

%% Constants init
cst = constants_init();

%% Parameters init (to default values, see values within routine)
default = parameters_default_init(cst);