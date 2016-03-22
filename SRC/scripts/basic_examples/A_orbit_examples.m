%--------------------------------------------------------------------------
% Example n°1: this matlab file makes use of the richardson third order
% approximation to build an EML1 planar lyapunov orbit (~8000 km of radius)
% Then, an EML2 halo orbit is produced (~10 000km of vertical extension)
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Orbit init & computation
lyap = init_lyap_orbit(cr3bp, cr3bp.l1, 8000 , cst);
lyap = lyapunov_orbit_computation(cr3bp, lyap, default, cst);

%% Same for a halo orbit
halo = init_halo_orbit(cr3bp, cr3bp.l2, cst.orbit.NORTHERN, 10000, cst);
halo = halo_orbit_computation(cr3bp, halo, default, cst);