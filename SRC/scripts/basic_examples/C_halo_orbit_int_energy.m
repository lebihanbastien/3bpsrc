%--------------------------------------------------------------------------
% Example n°3: 
%
% 1. Makes of the abacus
% ./data/halo_init_matrix_EML2.dat to generate an EML2 halo orbit 
% with a given energy value
%
% 2. Makes of the abacus
% ./data/halo_init_matrix_EML1.dat to generate an EML1 halo orbit 
% with the same energy value
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
%Orbit
orbit_L1 = init_halo_orbit_energy(cr3bp, cr3bp.l1, cst.orbit.NORTHERN, 3.12, cst);
orbit_L2 = init_halo_orbit_energy(cr3bp, cr3bp.l2, cst.orbit.NORTHERN, 3.12, cst);


%% Orbit computation
orbit_L1 = halo_orbit_interpolation_energy(cr3bp, orbit_L1, halo_init_EML1, default, cst);
orbit_L2 = halo_orbit_interpolation_energy(cr3bp, orbit_L2, halo_init_EML2, default, cst);
