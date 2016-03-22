%-------------------------------------------------------------------------%
% Computation of orbit.y0 and various orbital parameters from a third-order 
% guess (Richardson, 1980).
%
% WARNING 2: the use of third-order approximation should be limited to 
% Az < 30000km
%-------------------------------------------------------------------------%
function orbit = halo_orbit_computation(cr3bp, orbit, params, cst)

%-------------------------------------------------------------------------%
% Init from third order
%-------------------------------------------------------------------------% 
%Third order approximation
yv_to =  third_order_orbit(orbit, 0.0, cst);

% Integration vector
yv0 = (1:42)';

% 6-dim state
for k =1:6 
     yv0(k) = yv_to(k);   
end

% Specific case of the SOUTHERN family
if(strcmp(orbit.family,cst.orbit.family.SOUTHERN))
    yv0(3) = -yv0(3);    
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
if(orbit.Az_estimate > 0.0520 && params.diff_corr.type ~= cst.corr.X0_FIXED)
    params.diff_corr.type = cst.corr.X0_FIXED;
    disp('WARNING: the differential correction type has been changed because the desired Az is big (see params)');
end

orbit = diff_corr_3D(orbit.y0, cr3bp, orbit, params, cst);

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
    if(params.plot.halo_orbit == cst.TRUE)
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
if(params.plot.halo_orbit == cst.TRUE) %plotting  
  orbit_plot(yv, orbit, params, cst);
end

end