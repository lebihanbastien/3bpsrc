%-------------------------------------------------------------------------%
% Lyapunov orbits:
% Computation of orbit.y0 and various orbital parameters from a user-provided 
% guess 
% WARNING: the use of third-order approximation should be limited to
% Ax < 20000km
%-------------------------------------------------------------------------%
% @param cr3bp the structure containing the parent CR3BP
% @param orbit the structure containing the parent orbit
% @param params the structure containing the comptuation parameters
% @param cst the structure containing the numerical constants
% @param yv_to the user-provided guess of the initial conditions
% @return the udpated structure orbit
function orbit = lyapunov_orbit_continuation(cr3bp, orbit, params, cst, yv_to)

% Integration vector
yv0 = (1:42)';

% 6-dim state
for k =1:6 
     yv0(k) = yv_to(k);   
end

%STM concatenation after the 6-dim state
for i = 1 : 6
    for j = 1 : 6
        m = 6*(i-1) + j;
        yv0(m+6) = cst.orbit.STM0(i,j);
    end
end

%First guess (output if no diff correction)
orbit.y0 = yv0;

%-------------------------------------------------------------------------%
% Differential correction (if needed)
%-------------------------------------------------------------------------%
orbit = diff_corr_2D(orbit.y0, cr3bp, orbit, params);

%-------------------------------------------------------------------------%
% Period
%-------------------------------------------------------------------------%
% Integration over one 1/2 orbit
options = odeset('Events',@odezero_y,'Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
[~,~,te, yve] = ode45(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),[0 10],orbit.y0(1:6),options);
orbit.T12 = te; %1/2 period
orbit.T = 2*te; %period

%-------------------------------------------------------------------------%
% True Az
%-------------------------------------------------------------------------%
orbit.Az = max(abs(yve(3)), abs(yv0(3)));
orbit.Azdim = orbit.Az*cr3bp.L;

%-------------------------------------------------------------------------%
% Energy
%-------------------------------------------------------------------------%
orbit.C = jacobi(orbit.y0, cr3bp.mu);  %jacobi constant
orbit.E = -0.5*orbit.C;                %energy

%-------------------------------------------------------------------------%
% Linear algebra (monodromy matrix, etc)
%-------------------------------------------------------------------------%
%Integration over one orbit (42 variables: state + STM)
options = odeset('Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
[~,yv] = ode45(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),[0 orbit.T],orbit.y0,options);
    
%Monodromy matrix
orbit.monodromy = eye(6);
for i = 1 : 6
    for j = 1 : 6
        m = 6*(i-1) + j;
        orbit.monodromy(i,j) = yv(end,m+6);
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
if(params.plot.lyap_orbit == 1) %plotting 
    orbit_plot(yv, orbit, params, cst);
end
end