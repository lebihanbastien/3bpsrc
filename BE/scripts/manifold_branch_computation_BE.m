function manifold_branch = manifold_branch_computation_BE(cr3bp, orbit, manifold_branch, pos, t, params, cst)
% MANIFOLD_BRANCH_COMPUTATION computes a manifold leg associated to a given
% orbit in the CRTBP. Only UNSTABLE, INTERIOR manifolds are considered in
% this simplified version of the routine.
%
% MANIFOLD_BRANCH_COMPUTATION(CR3BP, ORBIT, MANIFOLD_BRANCH, POS, T,
% PARAMS, CST) computes a manifold leg in the system CR3BP associated to
% the orbit ORBIT, from information contained in the structure
% MANIFOLD_BRANCH. The starting position on the orbit is given by the
% variable POS (POS in [0 1] covers the entire orbit, but any value POS > 0
% is accepted). New initial conditions are taken along the unstable
% eigenvector, at a distance cr3bp.d_man from the orbit. The manifold leg 
% is then integrated up to the time T.
%
% See Koon et al. 2006, chapter 7, for details <a href="matlab:
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>.
%
%
% See also INIT_MANIFOLD_BRANCH
%
% BLB 2016

%--------------------------------------------------------------------------
% Integration on the orbit until the position POS is reached.
% POS should satisfy POS >= 0.
%--------------------------------------------------------------------------
if(pos > 0)
    [~, yv] = ode78_cr3bp([0 orbit.T*pos], orbit.y0, cr3bp.mu);
else
    yv = orbit.y0';
end

%--------------------------------------------------------------------------
% Current State Transition Matrix at POS
%--------------------------------------------------------------------------
%
% A COMPLETER (~ 1 Ligne)
%

%--------------------------------------------------------------------------
% New initial state at POS along the unstable direction
%--------------------------------------------------------------------------
%
% A COMPLETER (~ 5 Lignes)
%
% Rq1 : n'oubliez pas que vous avez déjà la direction instable à une
% position donnée de l'orbite grâce à orbit_postprocess.
% Rq2 : Le nouvel état le long de la direction instable doit s'appeler xs01

%--------------------------------------------------------------------------
% Save initial point in manifold
%--------------------------------------------------------------------------
manifold_branch.yv0 = xs01;

%--------------------------------------------------------------------------
% Integration time span
%--------------------------------------------------------------------------
tspan = [0 t];

%--------------------------------------------------------------------------
% Integration
%--------------------------------------------------------------------------
[te, yve, ~, ytraj] = ode78_cr3bp(tspan, xs01, cr3bp.mu);

%--------------------------------------------------------------------------
% Output
%--------------------------------------------------------------------------
manifold_branch.te    = te;      %time of the events
manifold_branch.yve   = yve;     %state at the events
manifold_branch.ytraj = ytraj;   %entire trajectory

%--------------------------------------------------------------------------
% Plotting (potentially)
%--------------------------------------------------------------------------
if(params.plot.manifold_branch)
    manifold_plot(ytraj, orbit, manifold_branch, params, cst);
end


end