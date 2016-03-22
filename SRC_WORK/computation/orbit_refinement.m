%-------------------------------------------------------------------------%
% orbit_refinement(cr3bp, orbit, params, cst)
%
% Computation of the initial conditions orbit.y0 and various orbital 
% parameters from a given guess yv_guess
%
% The routine makes use of a differential corrector adapted to three types
% of orbits: halo, planar lyapunov and vertical lyapunov. As for all
% differential correctors, the convergence strongly depends on the quality
% of the supplied first guess.
%-------------------------------------------------------------------------%
% Inputs:
% 1. cr3bp the structure containing the parent CR3BP
% 2. orbit the structure containing the orbit to update
% 3. params the structure containing the computation parameters
% 4. yv_guess the initial guess
% 5. cst the structure containing the numerical constants
%
% Outputs:
% 1. the udpated structure orbit
%
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function orbit = orbit_refinement(cr3bp, orbit, params, yv_guess, cst, isOrbitOnly)
               
% Integration vector
yv0 = (1:42)';

% 6-dim state
yv0(1:6) = yv_guess(1:6);   

%STM concatenation after the 6-dim state
yv0 = matrixToVector(yv0, cst.orbit.STM0, 6, 6, 6);

%First guess (output if no diff correction)
orbit.y0 = yv0;

%-------------------------------------------------------------------------%
% Change of differential correction procedure if necessary
% Base on heuristics on the Earth-Moon problem, not very robust!
%-------------------------------------------------------------------------%
if(isequal(orbit.type, cst.orbit.type.HALO) || isequal(orbit.type,cst.orbit.type.VLYAP))
    if(isfield(orbit, 'Az_estimate') == 1 && orbit.Az_estimate > 0.0520 && params.diff_corr.type ~= cst.corr.X0_FIXED)
        params.diff_corr.type = cst.corr.X0_FIXED;
        %disp('WARNING: the differential correction type has been changed because the desired Az is big (see params)');
    end
end

%-------------------------------------------------------------------------%
% Differential correction. Two steps:
% 
%   1. Select the dimensions of the differential correction procedure (dcp):
%       - The elements of the initial conditions (x0, y0...) that will be
%       modified. The corresponding indices are stored in xi0
%       - The elements of the final state (xf, yf, ...) that are targeted.
%       The corresponding indices are stored in xif.
%       - The section si = 0 (usually y = 0), on which the integration will
%       be stopped.
%
%   2. Perform the dcp.
%-------------------------------------------------------------------------%
%Switch between orbit types
switch(orbit.type)
    case cst.orbit.type.HALO
       
        %---------------------------------------------
        % Select the dimensions of the dcp
        %---------------------------------------------
        if(isequal(params.diff_corr.type, cst.corr.X0_FIXED))
            xi0 = [3, 5];  %(z0, vy0) are corrected
        elseif(isequal(params.diff_corr.type, cst.corr.Z0_FIXED))
            xi0 = [1, 5];  %(x0, vy0) are corrected
        else
            disp('orbit_refinement. Wrong type of differential correction for Halo orbits. return.');
        end
        xif = [4, 6]; %(vx0, vz0) = 0 is targeted
        si = 2; %the termination section is y = 0
        
        %---------------------------------------------
        % Perform the dcp
        %---------------------------------------------
        orbit = diff_corr_3D(orbit.y0, cr3bp , orbit, xi0, xif, si, params, cst);
        
        %---------------------------------------------
        % Old version
        %---------------------------------------------
        %orbit = diff_corr_HALO(orbit.y0, cr3bp, orbit, params, cst);
        
    case cst.orbit.type.VLYAP
        %---------------------------------------------
        % Select the dimensions of the dcp
        %---------------------------------------------
        if(isequal(params.diff_corr.type, cst.corr.X0_FIXED))
            xi0 = [3, 5];  %(z0, vy0) are corrected
        elseif(isequal(params.diff_corr.type, cst.corr.Z0_FIXED))
            xi0 = [1, 5];  %(x0, vy0) are corrected
        else
            disp('orbit_refinement. Wrong type of differential correction for Halo orbits. return.');
        end
        xif = [3, 4]; %(z0, vx0) = 0 is targeted
        si = 2;       %the termination section is y = 0
        
        %---------------------------------------------
        % Perform the dcp
        %---------------------------------------------
        orbit = diff_corr_3D(orbit.y0, cr3bp , orbit, xi0, xif, si, params, cst);
        
        %---------------------------------------------
        % Old version
        %---------------------------------------------
        %orbit = diff_corr_VLYAP(orbit.y0, cr3bp, orbit, params, cst);
        
        
    case cst.orbit.type.PLYAP
        orbit = diff_corr_2D(orbit.y0, cr3bp, orbit, params, cst);
end


%-------------------------------------------------------------------------%
% Period
%-------------------------------------------------------------------------%
if(params.computation.type == cst.computation.MATLAB)
    %-----------------------------
    % If MATLAB routines only
    %-----------------------------
    options = odeset('Events',@odezero_y,'Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
    [~,~,te, yve] = ode45(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),[0 10],orbit.y0(1:6),options);
else
    %-----------------------------
    % If MEX routines are allowed
    %-----------------------------
    %Event structure, for y = 0 section
    val_par = init_event(cst.manifold.event.type.Y_SECTION,...
        0.0,...
        cst.manifold.event.isterminal.YES,...
        cst.manifold.event.direction.ALL,...
        cr3bp.m1.pos,...
        cst);
    %Integration over one 1/2 orbit
    [te, yve] = ode78_cr3bp_event(0.0, 10, orbit.y0(1:6), 6, cr3bp.mu, val_par);
end

%Switch between orbit types for 1/2 period
switch(orbit.type)
    case {cst.orbit.type.HALO, cst.orbit.type.PLYAP}
        orbit.T12 = te(end);   %1/2 period
    case cst.orbit.type.VLYAP
        orbit.T12 = 2*te(end); %1/2 period
end
orbit.T = 2*orbit.T12; %period

%-------------------------------------------------------------------------%
% True Az/Ax
%-------------------------------------------------------------------------%
switch(orbit.type)
    case {cst.orbit.type.HALO, cst.orbit.type.VLYAP}
        orbit.Az    = max(abs(yve(3)), abs(yv0(3)));
        orbit.Azdim = orbit.Az*cr3bp.L;
    case cst.orbit.type.PLYAP
        orbit.Ax    = max(abs(yve(1)), abs(yv0(1)));
        orbit.Axdim = orbit.Ax*cr3bp.L;
end


%-------------------------------------------------------------------------%
% Energy
%-------------------------------------------------------------------------%
orbit.C = jacobi(orbit.y0, cr3bp.mu);  %jacobi constant
orbit.E = -0.5*orbit.C;                %energy

%-------------------------------------------------------------------------%
% Linear algebra (monodromy matrix, etc)
%-------------------------------------------------------------------------%
%Integration over one orbit (42 variables: state + STM)
if(params.computation.type == cst.computation.MATLAB)
    %-----------------------------
    % If MATLAB routines only
    %-----------------------------
    options = odeset('Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
    [~,yv] = ode45(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),[0 orbit.T],orbit.y0, options);
    yf = yv(end,:);
else
    %-----------------------------
    % If MEX routines are allowed
    %-----------------------------
    if(params.plot.orbit == cst.TRUE)
        [~, yf, ~, yv] = ode78_cr3bp(0.0, orbit.T, orbit.y0, 42, cr3bp.mu);
    else
        [~, yf] = ode78_cr3bp(0.0, orbit.T, orbit.y0, 42, cr3bp.mu);
    end
end

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

%Eigenvalures
for i = 1:6
    orbit.eigenvalues(i) = E(i,i);
end

%Stable and unstable direction (linear approx of the manifolds)
[~, posEigen] = min(orbit.eigenvalues);
orbit.stable_direction = V(:,posEigen);
[~, posEigen] = max(orbit.eigenvalues);
orbit.unstable_direction = V(:,posEigen);

%-------------------------------------------------------------------------%
% Status
%-------------------------------------------------------------------------%
orbit.status = cst.orbit.REAL;

%-------------------------------------------------------------------------%
% Plotting (potentially)
%-------------------------------------------------------------------------%
if(params.plot.orbit == 1) %plotting 
    orbit_plot(yv, orbit, params, cst, isOrbitOnly);
end
    
end