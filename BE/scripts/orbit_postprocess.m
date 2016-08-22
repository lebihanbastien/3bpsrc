function orbit = orbit_postprocess(cr3bp, orbit)
%   ORBIT_POSTPROCESS Postprocessing routine for generic orbit.
%
%   ORBIT_POSTPROCESS(CR3BP, ORBIT, PARAMS, CST) computes some key elements
%   of the orbit. Namely:
%
%   - Either the couple (Az, Azdim) - vertical extension for halo and
%   vertical orbits, or the couple (Ax, Axdim), maximum planar extension
%   for planar lyapunov orbits.
%   - orbit.C: the jacobian constant
%   - orbit.E: the energy
%   - orbit.yv: the state along the orbit on a given grid over the interval
%   [0, orbit.T].
%   - orbit.monodromy: the monodromy matrix.
%   - orbit.eigenvalues: the eigenvalues of the monodromy matrix in vector
%   form.
%   - orbit.stable_direction: the stable eigenvector of the monodromy
%   matrix.
%   - orbit.unstable_direction: the unstable eigenvector of the monodromy
%   matrix.
%
% See Koon et al. 2006, chapter 6, for details <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>.
%
%   This routine makes the following field of the orbit structure:
%   - orbit.y0:         initial conditions.
%   - orbit.T12:        half period.
%   - orbit.T:          period.
%
%   BLB 2016

%--------------------------------------------------------------------------
% True Az/Ax
%--------------------------------------------------------------------------
%Integration over one half orbit
[~, yve] = ode78_cr3bp([0 orbit.T12], orbit.y0(1:6), cr3bp.mu);


% Choose the maximum value between yve and yv0.
orbit.Az    = max(abs(yve(3)), abs(orbit.y0(3)));
orbit.Azdim = orbit.Az*cr3bp.L;


%--------------------------------------------------------------------------
% Energy
%--------------------------------------------------------------------------
orbit.C = jacobi(orbit.y0, cr3bp.mu);  %jacobi constant
orbit.E = -0.5*orbit.C;                %energy

%--------------------------------------------------------------------------
% Integration over one orbit (42 variables: state + STM)
%--------------------------------------------------------------------------
[~, yf, ~, yv] = ode78_cr3bp([0 orbit.T], orbit.y0, cr3bp.mu);


%--------------------------------------------------------------------------
% Save the state on a grid over one period:
%   - orbit.yv: the state along the orbit on a given grid over the interval
%   [0, orbit.T].
%--------------------------------------------------------------------------
orbit.yv = yv;

%--------------------------------------------------------------------------
% Linear algebra. Computation of: 
%   - orbit.monodromy: the monodromy matrix.
%   - orbit.eigenvalues: the eigenvalues of the monodromy matrix in vector
%   form.
%   - orbit.stable_direction: the stable eigenvector of the monodromy
%   matrix.
%   - orbit.unstable_direction: the unstable eigenvector of the monodromy
%   matrix.
%--------------------------------------------------------------------------
%
% A COMPLETER (~ 10-15 Lignes)
%

end
