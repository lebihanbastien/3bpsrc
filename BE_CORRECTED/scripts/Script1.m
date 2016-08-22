% First script of the BE.
% Computing a Halo orbit around EML2.
%
% The objective of this script is to build a periodic 3D halo orbit around
% the Libration point 2 of the Earth-Moon system (EML2).
%
% 
% Routines to complete: 
%      - diff_corr_3D_bb
%      - orbit_postprocess

%% Initialization of the workspace, constants and default parameters.
addpath(genpath('../'));
init;

%% Initialization of the environment:
% cr3bp is a matlab structure corresponding the the Earth-Moon Circular
% Restricted Three-Body Problem.
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Initialization of the orbit
orbit = init_orbit(cr3bp, ...                     % Parent CR3BP
                   cr3bp.l2, ...                  % Parent libration point is EML2
                   cst.orbit.type.HALO, ...       % Halo orbit
                   cst.orbit.family.NORTHERN, ... % Northern class
                   12000, ...                     % Of vertical extension Az ~ 12000km
                   cst);                          % Numerical constants

%% Computation of the orbit: 
% The routines to complete are used inside the routine orbit_computation.
orbit = orbit_computation(cr3bp, orbit, default, cst);