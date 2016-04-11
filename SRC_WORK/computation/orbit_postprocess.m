function orbit = orbit_postprocess(cr3bp, orbit, params, cst)
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
if(params.computation.type == cst.computation.MATLAB)
    %-----------------------------
    % If MATLAB routines only
    %-----------------------------
    options = odeset('Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);
    [~,yv] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu), [0 orbit.T12], orbit.y0(1:6), options);
    yve = yv(end,:);
else
    %-----------------------------
    % If MEX routines are allowed
    %-----------------------------
    [~, yve] = ode78_cr3bp(0.0, orbit.T12, orbit.y0(1:6), 6, cr3bp.mu);
end

% Choose the maximum value between yve and yv0.
switch(orbit.type)
    case {cst.orbit.type.HALO, cst.orbit.type.VLYAP}
        orbit.Az    = max(abs(yve(3)), abs(orbit.y0(3)));
        orbit.Azdim = orbit.Az*cr3bp.L;
    case cst.orbit.type.PLYAP
        orbit.Ax    = max(abs(yve(1)), abs(orbit.y0(1)));
        orbit.Axdim = orbit.Ax*cr3bp.L;
end

%--------------------------------------------------------------------------
% Energy
%--------------------------------------------------------------------------
orbit.C = jacobi(orbit.y0, cr3bp.mu);  %jacobi constant
orbit.E = -0.5*orbit.C;                %energy

%--------------------------------------------------------------------------
% Change the period for Vertical orbits
%--------------------------------------------------------------------------
if(strcmp(orbit.type, cst.orbit.type.VLYAP))
    orbit.T12 = 2*orbit.T12;
    orbit.T   = 2*orbit.T;
end

%--------------------------------------------------------------------------
%Integration over one orbit (42 variables: state + STM)
%--------------------------------------------------------------------------
if(params.computation.type == cst.computation.MATLAB)
    %-----------------------------
    % If MATLAB routines only
    %-----------------------------
    options = odeset('Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);
    [~, yv] = ode113(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),[0 orbit.T],orbit.y0, options);
    yf = yv(end,:);
else
    %-----------------------------
    % If MEX routines are allowed
    %-----------------------------
    [~, yf, ~, yv] = ode78_cr3bp(0.0, orbit.T, orbit.y0, 42, cr3bp.mu);
end

%--------------------------------------------------------------------------
% Save the state on a grid over one period
%--------------------------------------------------------------------------
orbit.yv = yv;

%--------------------------------------------------------------------------
% Max/Min distance with each primary
%--------------------------------------------------------------------------
[orbit.minDistToM1, ~, orbit.maxDistToM1]  = distToPrimary(yv, cr3bp.m1);
[orbit.minDistToM2, minDistIndex, orbit.maxDistToM2, maxDistIndex] = distToPrimary(yv, cr3bp.m2);

%--------------------------------------------------------------------------
% 'Perigee' wrt M2 (Moon in the Earth-Moon system)
%--------------------------------------------------------------------------
orbit.perigee.radius   = orbit.minDistToM2;
orbit.perigee.altitude = orbit.minDistToM2-cr3bp.m2.Rm/cr3bp.L;
orbit.perigee.position = yv(minDistIndex, 1:3);

%--------------------------------------------------------------------------
% 'Apogee' wrt M2 (Moon in the Earth-Moon system)
%--------------------------------------------------------------------------
orbit.apogee.altitude = orbit.maxDistToM2;
orbit.apogee.position = yv(maxDistIndex, 1:3);

%--------------------------------------------------------------------------
% Linear algebra (monodromy matrix, etc)
%--------------------------------------------------------------------------
%Monodromy matrix
orbit.monodromy = eye(6);
for i = 1 : 6
    for j = 1 : 6
        m = 6*(i-1) + j;
        orbit.monodromy(i,j) = yf(m+6);
    end
end

%Eigen
[V,E] = eig(orbit.monodromy);

%Eigenvalues
for i = 1:6
    orbit.eigenvalues(i) = E(i,i);
end

%Stable and unstable direction (linear approx of the manifolds)
[~, posEigen] = min(orbit.eigenvalues);
orbit.stable_direction = V(:,posEigen);
[~, posEigen] = max(orbit.eigenvalues);
orbit.unstable_direction = V(:,posEigen);

end


function [minDist, minDistIndex, maxDist, maxDistIndex] = distToPrimary(yv, primary)
% [MINDIST, MINDISTINDEX, MAXDIST, MAXDISTINDEX] = DISTTOPRIMARY(YV, PRIMARY)
% computes the the min/max distance between the current state
% YV and the center of the primary PRIMARY.
%
% BLB 2016

% Get the difference in position for each values in yv
ydist  = bsxfun(@minus,yv(:,1:3), primary.pos);
% Get the norm of each difference vector
yrdist = arrayfun(@(idx) norm(ydist(idx,:)), 1:size(ydist,1));
% Get the mininum
[minDist, minDistIndex] = min(yrdist);
% Get the maximum
[maxDist, maxDistIndex] = max(yrdist);
end