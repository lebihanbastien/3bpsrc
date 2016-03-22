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
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
%Orbit
orbit = init_halo_orbit(cr3bp, cr3bp.l2, cst.orbit.NORTHERN, 30000, cst);
%Interpolation matrix
halo_init = halo_init_EML2;

%% Orbit computation
orbit = halo_orbit_interpolation(cr3bp, orbit, halo_init, default, cst);


%% Manifold computation
%Unstable manifold
%-------------------------
manifold_branch_unstable_exterior  = init_manifold_branch_event(cst.manifold.UNSTABLE,...
                                          cst.manifold.EXTERIOR,...
                                          cst.manifold.event.type.X_SECTION,...
                                          0.0,...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,... 
                                          cr3bp.m1.pos,...
                                          cst);
                                      
manifold_branch_unstable_interior  = init_manifold_branch_event(cst.manifold.UNSTABLE,...
                                          cst.manifold.INTERIOR,...
                                          cst.manifold.event.type.X_SECTION,...
                                          0.0,...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,...
                                          cr3bp.m1.pos,...
                                          cst);

%Stable manifold
%-------------------------
manifold_branch_stable_exterior  = init_manifold_branch_event(cst.manifold.STABLE,...
                                          cst.manifold.EXTERIOR,...
                                          cst.manifold.event.type.X_SECTION,...
                                          0.0,...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,...
                                          cr3bp.m1.pos,...
                                          cst);
                                      
manifold_branch_stable_interior  = init_manifold_branch_event(cst.manifold.STABLE,...
                                          cst.manifold.INTERIOR,...
                                          cst.manifold.event.type.X_SECTION,...
                                          0.0,...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,...
                                          cr3bp.m1.pos,...
                                          cst);
                                      
                                      
%% Computation and plot
% Integration duration
t = 10;


% Unstable interior
for theta = 0:0.05:1 %position on the orbit
    manifold_branch_unstable_interior = manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable_interior, theta, t, default, cst);
end

% Unstable exterior
for theta = 0:0.05:1 %position on the orbit
    manifold_branch_unstable_exterior = manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable_exterior, theta, t, default, cst);end

% Stable exterior
for theta = 0:0.05:1 %position on the orbit
    manifold_branch_stable_exterior = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable_exterior, theta, t, default, cst);end

% Stable interior
for theta = 0:0.05:1 %position on the orbit
    manifold_branch_stable_interior = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable_interior, theta, t, default, cst);
end


%% Print
orbit = halo_orbit_interpolation(cr3bp, orbit, halo_init, default, cst);
print('manifold_stable_interior.eps', figure(1), '-depsc')
print('manifold_stable_interior_3d.eps', figure(4), '-depsc')

print('manifold_stable_interior.png', figure(1), '-dpng')
print('manifold_stable_interior_3d.png', figure(4), '-dpng')


