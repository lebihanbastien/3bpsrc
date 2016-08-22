%--------------------------------------------------------------------------
% This matlab file makes of the abacus
% ./data/halo_init_matrix_EML1.dat to generate an EML1 halo orbit and its
% interior stable and unstable manifolds.
%
% WARNING 1: The abacuses provided here are only valid for the Earth-Moon
% Lagrange points 1 & 2. For other systems, one needs to use the routine 
% halo_orbit_computation, instead of halo_orbit_interpolation.
%
% WARNING 2: a good plot of the manifolds greatly depends on an adapted
% integration time wrt the orbit size (i.e. a big orbit calls for a big
% integration time to be able to see a divergence from the original orbit)
%
% BLB 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Change default parameter
default.plot.firstPrimDisp = false;

%% Structures init
%Environment
cr3bp = init_CR3BP('SUN', 'EARTH', default);
%Orbit
orbit = init_orbit(cr3bp, ...                     % Parent CR3BP
                   cr3bp.l1, ...                  % Parent libration point
                   cst.orbit.type.HALO, ...       % HALO orbit 
                   cst.orbit.family.SOUTHERN, ... % Southern class
                   400000, ...                     % Of vertical extension ~ 30000 km
                   cst);                          % Numerical constants
               

%% Orbit computation
orbit = orbit_computation(cr3bp, orbit, default, cst);

%% Manifold initialization

% We define a stable manifold
manifold_branch_stable    = init_manifold_branch(cst.manifold.STABLE, ...
                                                 cst.manifold.EXTERIOR);
% We define an unstable manifold                                            
manifold_branch_unstable  = init_manifold_branch(cst.manifold.UNSTABLE,...
                                                 cst.manifold.EXTERIOR);
                                             
% We define an event structure that can trigger the end of the integration 
% along the manifold. Here, the integration will stop when the line that
% links the current state and the center of the Earth reaches a given angle
% with respect to the Earth-Moon line.
earth.event = init_event(cst.manifold.event.type.X_SECTION,...         %the event is triggered at a given angle...
                        cr3bp.m2.pos(1),...                                    %of a given value...
                        cst.manifold.event.isterminal.YES,...              %the trajectory stops at the first ocurrence...
                        cst.manifold.event.direction.INCREASING,...               %all direction are considered...
                        cr3bp.m2.pos, cst);                                %the center for the computation of the angle is the Earth

% Note: the event can be directly incorporated to the manifold during its
% initialization (see init_manifold_branch.m).

%% Manifold computation                                             
% Integration duration
t = 10;

% Stable
for theta = 0:0.05:1 %position on the orbit
    manifold_branch_stable = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable, theta, t, default, cst, earth.event);
end

% Unstable, with an event triggered.
for theta = 0:0.05:1 %position on the orbit
    manifold_branch_unstable = manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable, theta, t, default, cst, earth.event);
end

% Note: the routine manifold_branch_computation allows the event to be
% given as a structure or as a function handler (see
% manifold_branch_computation.m). In this case, the equivalent call would
% be: 
% manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable,...
%                             theta, t, default, cst,...
%                             @(t,y)odezero_angle(t,y,earth.event))
% The function handlers allow the user to define its own event
% routines, on top of the usual ones, defined in the structure
% cst.manifold.event.type

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end