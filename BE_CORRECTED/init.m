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

%% Constants init
cst = constants_init();

%% Parameters init (to default values, see values within routine)
default = parameters_default_init(cst);