%--------------------------------------------------------------------------
% Continuation example n°2: Continuation procedure to produce a discrete set within 
% the family of northern orbits of EML2
% WARNING: very sensitive, works only for small Az amplitudes (contrary to
% planar lyapunov case). Complete continuations are more complex than that
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Orbit
orbit_1 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.HALO, cst.orbit.family.NORTHERN, 8000, cst);
orbit_2 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.HALO, cst.orbit.family.NORTHERN, 10000, cst);

%% Orbit computation
orbit_1 = orbit_computation(cr3bp, orbit_1, default, cst);
orbit_2 = orbit_computation(cr3bp, orbit_2, default, cst);

%% Continuation
% CAREFUL: very sensitive, works only for small Az amplitudes (contrary to
% planar lyapunov case). Complete continuations are more complex than that
for i = 1:10
    yv = orbit_2.y0 + (orbit_2.y0 - orbit_1.y0);
    orbit_1 = orbit_2;
    orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yv, cst);
end

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end


