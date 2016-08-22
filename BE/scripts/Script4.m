% Fourth script of the BE.
% Computing a Halo orbit around EML2 + its unstable directions + the
% complete unstable manifold, stopping at the closest approach of the lunar
% surface.
%
% Objectives: 
%  - Produce the tube-like shape of the unstable manifold.
%  - See the influence of the size of the orbit (Az) on the reachable
% altitudes with respect to the Lunar surface. Az must be taken below
% 20000 km
%
% 2016

%% Call of the first script: building the orbit
Script1

%% Stop plotting on figure 1
default.plot.XY = false;

%% Integration duration of the manifold arbitrarily fixed to 22 days
t = 22*cst.env.days*2*pi/cr3bp.T;

%% Initialization of the manifold
msi = init_manifold_branch(cst.manifold.UNSTABLE, cst.manifold.INTERIOR);

%% Initialize some useful constants
% Lunar radius
rm = cr3bp.m2.Rm/cr3bp.L;
% Lunar position
pm = cr3bp.m2.pos;
% Earth-Moon distance
Lf = cr3bp.L;

%% Manifold computations
%
% A COMPLETER (~ 10-15 Lignes)
%

%% Plot
%
% A COMPLETER (~ 5 Lignes)
%