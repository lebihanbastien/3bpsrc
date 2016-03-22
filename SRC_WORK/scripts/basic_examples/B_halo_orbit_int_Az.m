%--------------------------------------------------------------------------
% Example n°2: this matlab file makes of the abacus
% ./data/halo_init_matrix_EML1.dat to generate an EML1 halo orbit and its
% interior stable and unstable manifolds.
%
% WARNING 1: The abacuses provided here are only valid for the Earth-Moon
% Lagrange points 1 & 2. For other systems, one needs to use the routine 
% halo_orbit_computation, instead of halo_orbit_interpolation.
%
% WARNING 2: a good plot of the manifolds greatly depends on an adapted
% integration time wrt the orbit size (i.e. a big orbit calls for a big
% integration time to be able to see a divergence from the original motion)
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------

% Reboot
clear all;
close all;

% Add subfolders to the path
addpath('./computation');
addpath('./data');
addpath('./init');
addpath('./ode');

%% Data loading (abacus)
load halo_init_matrix_EML2 halo_init_EML2;
load halo_init_matrix_EML1 halo_init_EML1;

%% Constants init
cst = constants_init();

%% Parameters init (to default values, see values within routine)
default = parameters_default_init(cst);

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
%Orbit
orbit = init_orbit(cr3bp, ...                     % Parent CR3BP
                   cr3bp.l1, ...                  % Parent libration point
                   cst.orbit.type.HALO, ...       % HALO orbit 
                   cst.orbit.family.SOUTHERN, ... % Southern class
                   30000, ...                     % Of vertical extension ~ 30000 km
                   cst);                          % Numerical constants
               
%Interpolation matrix
halo_init = halo_init_EML1;

%% Orbit computation
orbit = halo_orbit_interpolation(cr3bp, orbit, halo_init, default, cst);

%% Manifold initialization
manifold_branch_stable    = init_manifold_branch(cst.manifold.STABLE, ...
                                                 cst.manifold.INTERIOR, cst);
                                             
manifold_branch_unstable  = init_manifold_branch(cst.manifold.UNSTABLE,...
                                                 cst.manifold.INTERIOR, cst);

%% Manifold computation                                             
% Integration duration
t = 5;

% Stable
for theta = 0:0.1:1 %position on the orbit
    manifold_branch_stable = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable, theta, t, default, cst);
end

% Unstable
for theta = 0:0.1:1 %position on the orbit
    manifold_branch_unstable = manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable, theta, t, default, cst);
end

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end