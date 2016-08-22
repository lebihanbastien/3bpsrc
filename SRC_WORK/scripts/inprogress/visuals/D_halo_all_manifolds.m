%--------------------------------------------------------------------------
% Example n°4: this matlab file makes of the abacus
% ./data/halo_init_matrix_EML2.dat to generate an EML2 halo orbit and its
% exterior stable and unstable manifolds starting with a given Az
%
% Then, the stable/unstable interior/exterior manifolds are computed up to 
% the section x = 0. Note that, because of the gravitationnal influence of
% the Moon, some legs of the interior manifolds can appear to be on the
% corresponding exterior manifolds, and vice versa.
%
% WARNING: a good plot of the manifolds greatly depends on an adapted
% integration time to the orbit size (i.e. a big orbit calls for a big
% integration time to be able to see a divergence from the original motion)
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters

default.plot.XY              = false; %plot also the results in X-Z plane
default.plot.firstPrimDisp   = cst.TRUE;  %is the first primary (e.g. the Sun in the Sun-Earth system) displayed?
default.plot.allLibPoints    = cst.TRUE;  %are all libration points displayed?
default.plot.names           = cst.TRUE;  %are the names displayed?
default.plot.tdAxes          = cst.TRUE;  %are the pretty 3D axes displayed?
default.plot.bigPrimFac      = 3.0;       %the primaries appear bigPrimFac x bigger than they actually are (easier to see on screen)

% 4. See parameters_default_init.m to see other options
%--------------------------------------------------------------------------

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
%Orbit
orbit = init_orbit(cr3bp, ...                     % Parent CR3BP
                   cr3bp.l1, ...                  % Parent libration point
                   cst.orbit.type.HALO, ...       % HALO orbit 
                   cst.orbit.family.NORTHERN, ... % Northern class
                   10000, ...                     % Of vertical extension ~ 30000 km
                   cst);                          % Numerical constants
              
%Interpolation matrix
halo_init = halo_init_EML1;

%% Orbit computation
orbit = halo_orbit_interpolation(cr3bp, orbit, halo_init, default, cst);


%% Manifold computation
%Unstable manifold
%-------------------------
manifold_branch_unstable_exterior  = init_manifold_branch_event(cst.manifold.UNSTABLE,...
                                          cst.manifold.EXTERIOR,...
                                          cst.manifold.event.type.FREE,...
                                          degtorad(-60),...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,... 
                                          cr3bp.m1.pos,...
                                          cst);
                                      
manifold_branch_unstable_interior  = init_manifold_branch_event(cst.manifold.UNSTABLE,...
                                          cst.manifold.INTERIOR,...
                                          cst.manifold.event.type.FREE,...
                                          1-cr3bp.mu,...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,...
                                          cr3bp.m1.pos,...
                                          cst);

%Stable manifold
%-------------------------
manifold_branch_stable_exterior  = init_manifold_branch_event(cst.manifold.STABLE,...
                                          cst.manifold.EXTERIOR,...
                                          cst.manifold.event.type.FREE,...
                                          degtorad(60),...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,...
                                          cr3bp.m1.pos,...
                                          cst);
                                      
manifold_branch_stable_interior  = init_manifold_branch_event(cst.manifold.STABLE,...
                                          cst.manifold.INTERIOR,...
                                          cst.manifold.event.type.FREE,...
                                          1-cr3bp.mu,...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,...
                                          cr3bp.m1.pos,...
                                          cst);
                                      
                                      
%% Computation and plot
% Integration duration
t = 10;

%%
%Unstable interior
for theta = 0:0.1:1 %position on the orbit
    manifold_branch_unstable_interior = manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable_interior, theta, t, default, cst);
end

% Unstable exterior
for theta = 0:0.1:1 %position on the orbit
    manifold_branch_unstable_exterior = manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable_exterior, theta, t, default, cst);
end

% % Stable exterior
% for theta = 0:0.05:1 %position on the orbit
%     manifold_branch_stable_exterior = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable_exterior, theta, t, default, cst);
% end
% 
% % Stable interior
% for theta = 0:0.05:1 %position on the orbit
%     manifold_branch_stable_interior = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable_interior, theta, t, default, cst);
% end

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end

%% Go on with matlab 2015a
h = figure(4);
set(h, 'Color',[1 1 1]);
set(h, 'Position', [100, 100, 1050, 900]);
zoom(2)
hold on;
grid off;
axis off;