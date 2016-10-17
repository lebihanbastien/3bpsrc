%--------------------------------------------------------------------------
% Init script for all computation. Should be called at the beginning of
% each computation.
%
% BLB 2016
%--------------------------------------------------------------------------

%% Reboot
clear all;
close all;

%% Plot settings
set(0, 'defaultTextFontSize', 18);
set(0, 'defaultAxesFontSize', 18);
set(0, 'defaultTextFontWeight', 'bold');
set(0, 'defaultTextHorizontalAlignment', 'center');
set(0, 'defaultLineMarkerSize', 2);

%% Add subfolders to the path
% Recursively add all subfolders
addpath(genpath('.'));

%% Data loading (abacuses)
load halo_init_matrix_EML2 halo_init_EML2;
load halo_init_matrix_EML1 halo_init_EML1;
% New abacuses (more complete)
load halo2_init_EML2 halo2_init_EML2;
load halo2_init_EML1 halo2_init_EML1;

load nro_init_EML1 nro_init_EML1;
load nro_init_EML2 nro_init_EML2;

%% Constants init
cst = constants_init();

%% Parameters init (to default values, see values within routine)
default = parameters_default_init(cst);