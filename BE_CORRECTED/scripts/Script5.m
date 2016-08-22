% Fifth script of the BE.
% Computing a Halo orbit around EML2 + its unstable directions + the
% complete unstable manifold, stopping at the closest approach of the lunar
% surface. Then, a study of the Low Lunar Orbits (LLOs) that are directly
% reachable with the unstable manifold is performed.
%
% At the closest approach to the lunar surface, we suppose that the
% spacecraft performs a maneuver tangent to its current velocity in order
% to set itself on a circular orbit around the Moon. This orbit has an
% altitude equal to the current altitude of the spacecraft. This altitude
% can be very large, and we are only interested in Low Lunar Orbits (LLOs)
% with an altitude below 500 km. 
%
% The objective of this script is to find the position on the initial halo
% orbit around EML2 that allows to reach a circular Low Lunar Orbit using
% the corresponding unstable manifold leg as the transfer trajectory.
%
% Ideally, we would also like to end up on "frozen" orbits around the Moon,
% which are orbits particularly stable. These orbits are found for
% an inclination i equal to 27°, 50°, 76°, and 86°.
%
% Objectives: 
%      - See how the altitude and inclination of the final lunar orbit are
%      related to the size of the initial halo orbit (Az) and the departure
%      position on this halo orbit.  Az must be taken below 20000 km.
%      - Find the Az that allows to get the first frozen orbits (i = 27°+-2°).
%

%% Call of the first script: building the orbit
Script1

%% Stop plotting on figure 1
default.plot.XY = false;

%% Initialization of the manifold
msi = init_manifold_branch(cst.manifold.UNSTABLE, cst.manifold.INTERIOR);
default.plot.manifold_branch = false;

%% Reachability of the Moon
% Vector of positions on the halo orbit.
thetav = 0:0.005:1;
% Maximum altitude allowed at the Moon : 500 km.
maxalt = 500/cr3bp.L; %careful: normalization!
% Moon reachability
[output] = moonreachability(cr3bp, orbit, msi, thetav, maxalt, default, cst);
