%--------------------------------------------------------------------------
% Example n°1: this matlab file makes use of the richardson third order
% approximation to build:
% - an EML2 planar lyapunov orbit (~8000 km of radius)
% - an EML2 halo orbit (~10 000km of vertical extension)
% - an EML2 vertical lyapunov orbit (~30 000km of vertical extension)
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters


%% Environment init
cr3bp = init_CR3BP('SUN', 'EARTH', default);


%% Same for a halo orbit
%Initialization
halo = init_orbit(cr3bp, ...       % Parent CR3BP
    cr3bp.l1, ...                  % Parent libration point
    cst.orbit.type.HALO, ...       % HALO orbit
    cst.orbit.family.NORTHERN, ... % Northern class
    177000, ...                      % Of vertical extension ~ 10000 km
    cst);                          % Numerical constants

%Computation
halo = orbit_computation(cr3bp, halo, default, cst);

%% Manifold initialization                                                                       
% The unstable manifold is stopped when unstable.event occurs
manifold_branch_unstable = init_manifold_branch(cst.manifold.UNSTABLE, cst.manifold.EXTERIOR);

%% Manifold computation
t = 336*86400*2*pi/cr3bp.T;

%% Stable
for theta = 0:0.1:0 %position on the orbit
    manifold_branch_unstable = manifold_branch_computation(cr3bp, halo, manifold_branch_unstable, theta, t, default, cst);
end
