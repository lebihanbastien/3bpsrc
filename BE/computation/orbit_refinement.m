function orbit = orbit_refinement(cr3bp, orbit, params, yvg, cst, varargin)
% ORBIT_REFINEMENT computation of halo
% symmetric periodic orbits in the CRTBP, from an initial guess.
%
% ORBIT = ORBIT_REFINEMENT(CR3BP, ORBIT, PARAMS, YVG, CST) computes an 
% orbit in the system CR3BP, with the desired characteristics
% listed in the structure ORBIT (size or energy, lagrange point), and from
% an first guess of initial condition YVG.
% A differential correction process is applied to get a real periodic 
% orbit.
% Finally, a post process routine is applied in order to compute various
% characteristics of the orbit, listed below. Depending on the user-defined
% parameter structure PARAMS, the result can be plotted at the end of the
% routine.
%
% See Koon et al. 2006, chapter 6, for details <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>.
%
% BLB 2016

%--------------------------------------------------------------------------
% Initialisation of the initial conditions:
% yv0 = [yvg ; State Transition Matrix]
%--------------------------------------------------------------------------
% Integration vector
yv0 = (1:42)';
% 6-dim state from a first approximation
yv0(1:6) = yvg(1:6);
% STM concatenation after the 6-dim state
yv0 = matrixToVector(yv0, cst.orbit.STM0, 6, 6, 6);

%--------------------------------------------------------------------------
% Differential correction procedure (DCP). At the end of this procedure,
% the following elements (at least) are updated in the orbit structure:
%   - orbit.y0:         initial conditions.
%   - orbit.T12:        half period.
%   - orbit.T:          period.
%--------------------------------------------------------------------------
orbit = diff_corr_3D_bb(yv0, cr3bp , orbit, params, cst);

%--------------------------------------------------------------------------
% Status
%--------------------------------------------------------------------------
orbit.status = cst.orbit.REAL;

end