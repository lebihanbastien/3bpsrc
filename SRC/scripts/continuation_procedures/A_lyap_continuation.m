%--------------------------------------------------------------------------
% Continuation example n°1: Continuation procedure to produce a discrete 
% set within the family of planar lyapunov orbits of EML2
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Orbit
orbit_1 = init_lyap_orbit(cr3bp, cr3bp.l2, 1000 , cst);
orbit_2 = init_lyap_orbit(cr3bp, cr3bp.l2, 2000 , cst);

%% Orbit computation
orbit_1 = lyapunov_orbit_computation(cr3bp, orbit_1, default, cst);
orbit_2 = lyapunov_orbit_computation(cr3bp, orbit_2, default, cst);

%% Continuation

for i = 1:50
    yv = orbit_2.y0 + (orbit_2.y0 - orbit_1.y0);
    orbit_1 = orbit_2;
    orbit_2 = lyapunov_orbit_continuation(cr3bp, orbit_2, default, cst, yv);
end



