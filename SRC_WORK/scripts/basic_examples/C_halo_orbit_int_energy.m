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
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters
%default.computation.type = cst.computation.MATLAB;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
%Orbit
orbit_L1 = init_halo_orbit_energy(cr3bp, cr3bp.l1, cst.orbit.family.NORTHERN, 3.12, cst);
orbit_L2 = init_halo_orbit_energy(cr3bp, cr3bp.l2, cst.orbit.family.NORTHERN, 3.12, cst);


%% Orbit computation
orbit_L1 = halo_orbit_interpolation(cr3bp, orbit_L1, halo_init_EML1, default, cst);
orbit_L2 = halo_orbit_interpolation(cr3bp, orbit_L2, halo_init_EML2, default, cst);

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end