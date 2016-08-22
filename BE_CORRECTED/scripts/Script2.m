% Second script of the BE.
% Computing a Halo orbit around EML2 + its unstable directions + a branch
% of the unstable manifold.
%
% The objective of this script is to compute the unstable direction on any
% position on the halo orbit of Script1. Then, initial conditions on this
% unstable direction are computed forward in time to form a branch (or leg)
% of the unstable manifold associated to this orbit.
%
% Routines to complete: 
%      - manifold_branch_computation_BE
%

%% Call of the first script to compute the orbit again
Script1

%% Stop plotting on Figure 1
default.plot.XY = false;

%% Integration duration: arbitrarily fixed to 22 days
t = 22*cst.env.days*2*pi/cr3bp.T;

%% Initialization of the manifold
msi = init_manifold_branch(cst.manifold.UNSTABLE,...    %it is an UNSTABLE manifold.
                           cst.manifold.INTERIOR);      %it is an INTERIOR manifold: it will go towards the Moon.
                       
%% Computation of the manifold
% Departure position on the orbit in the interval [0,1]
theta = 0.0;
% Computation of the manifold branch
msi = manifold_branch_computation_BE(cr3bp,...   %parent CRTBP
                                     orbit,...   %parent orbit
                                     msi,...     %current manifold branch
                                     theta,...   %departure position on the orbit
                                     t,...       %time of flight on the manifold branch
                                     default,... %default parameters
                                     cst);       %constants