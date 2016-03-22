%--------------------------------------------------------------------------
% Example n°6: Comparison of two different ways to compute halo orbits:
%
%   - The interpolation from a given abacus (halo_orbit_interpolation).
%   - The computation from a third-order first guess (halo_orbit_computation).
%
% This matlab file does two things:
%
%  1. It makes us of the abacus
% ./data/halo_init_matrix_EML1.dat and ./data/halo_init_matrix_EML2.dat
% to generate two halo orbits around L1 and L2 of same Az (15000km)
% via the routine:
%                   halo_orbit_interpolation
%  2. It computes the orbit of same Az_estimate starting from the 3rd
% order approximation of Richardson, via the routine:
%                   halo_orbit_computation
%--------------------------------------------------------------------------
% Important remarks:
%
% A. When the precision of the integration is high (see default.ode45
% structure), the difference in CPU time between the interpolation process 
% and the computation from an initial 3rd order guess can be important.
%
% B. For the Earth-Moon system, one can try to set Az = 35000 km. 
% For this value, the third order approximation of Richardson is clearly
% not valid anymore: the difference between the interpolated orbits (the
% good ones) and the orbits computed from the 3rd order approx is huge.
% For even higher values of Az, the computation based on Richardson's 3rd
% order diverges.
%
% C. However, the interpolation process requires abacuses such as 
% ./data/halo_init_matrix_EML1.dat in order to generate the solutions.
% in these version of the code, only the Earth-Moon L1,2 abacuses are
% provided. For other computation (either different points or different
% systems), the halo_orbit_computation routine MUST be used, keeping in
% mind that it is valid for only a certain range of Az values.
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% User data
user.Az = 15000; %Az in [km]

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Orbits
halo_L2 = init_halo_orbit(cr3bp, cr3bp.l2, cst.orbit.NORTHERN, user.Az, cst);
halo_L1 = init_halo_orbit(cr3bp, cr3bp.l1, cst.orbit.NORTHERN, user.Az, cst);

%% Orbit computation with interpolation
tic;
halo_L2 = halo_orbit_interpolation(cr3bp, halo_L2, halo_init_EML2, default, cst);
halo_L1 = halo_orbit_interpolation(cr3bp, halo_L1, halo_init_EML1, default, cst);
toc;

%% Orbit computation with third order approximation
tic;
halo_L2 = halo_orbit_computation(cr3bp, halo_L2, default, cst);
halo_L1 = halo_orbit_computation(cr3bp, halo_L1, default, cst);
toc;

