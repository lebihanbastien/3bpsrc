% Third script of the BE.
% Computing a Halo orbit around EML2 + its unstable directions + a branch
% of the unstable manifold, stopping at the closest approach of the lunar
% surface.
%
% There is nothing to do here. Just check that this script does basically
% the same result as Script2, except that now the integration of the
% manifold leg stops at the closest approach of the lunar surface.


%% Call of the first script.
Script1

%% Stop plotting on figure 1
default.plot.XY = false;

%% Integration duration: arbitrarily fixed to 22 days
t = 22*cst.env.days*2*pi/cr3bp.T;

%% Initialization of the manifold
msi = init_manifold_branch(cst.manifold.UNSTABLE, cst.manifold.INTERIOR);

%% Computation of a manifold that stops at the closest approach of the lunar surface
% Departure position on the orbit in [0, 1]
theta = 0.0;
% Computation of the manifold branch. The routine manifold_branch_computation_moon
% is equivalent to the routine manifold_branch_computation_BE but stops the
% integration at the closest approach of the lunar surface.
msi = manifold_branch_computation_moon(cr3bp, orbit, msi, theta, t, default, cst);