%--------------------------------------------------------------------------
% Plots an NRO orbit.
%
% BLB 2016
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Change of parameters wrt default
default.plot.XY          = true;  %plot also the results in X-Y plane
default.plot.XZ          = true;  %plot also the results in X-Z plane
default.plot.YZ          = true;  %plot also the results in Y-Z plane

%% Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Orbit
%Initialize NRO
nro = init_nro(cr3bp, cr3bp.l2, cst.orbit.family.SOUTHERN, cst);

% Interpolation and plot
nro = nro_interpolation(cr3bp, nro, nro_init_EML2, default, cst, 'altitudeOfPerigee', 1000);