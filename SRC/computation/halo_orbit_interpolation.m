function orbit = halo_orbit_interpolation(cr3bp, orbit, halo_init, params, cst)
%-------------------------------------------------------------------------%
% Polynomial interpolation in abacus to recover good guess for orbit.y0.
% The temporary structure fit contains all the elements for this interp.
%
% WARNING 1: the interpolation alone should give a good enough guess to
% directly produce a periodic orbit. The forcing of the differential 
% correction scheme is at the discretion of the user.
%
% WARNING 2: the interpolation is limited to: Az = 76698 km  for EML2
%                                             Az = 126788 km for EML1
%
% RM: the data files contain the NORTHERN family.
%-------------------------------------------------------------------------%
if(orbit.Az < halo_init.Azlimit)  
    fit.half = 2;
    fit.degree = 2*fit.half;
    [~, array_position] = min(abs(halo_init.matrix(:,7) - orbit.Az));  
    fit.x =  halo_init.matrix(array_position - fit.half: array_position + fit.half,7);
    for count =1:6
        fit.y(:,count) =  halo_init.matrix(array_position - fit.half: array_position + fit.half,count);
        %Fitting for every dimension of the state (6)
        [fit.p, ~, fit.mu] = polyfit(fit.x,fit.y(:,count),fit.degree);
        %Evaluation
        fit.f(count) = polyval(fit.p,orbit.Az,[],fit.mu);
    end
else
    disp('WARNING: the vertical extension Az is out of range in halo_init.matrix');
end

% Integration vector
yv0 = (1:42)';

% 6-dim state
for k =1:6 
     yv0(k) = fit.f(k);   
end

% Specific case of the SOUTHERN family
if(strcmp(orbit.family,cst.orbit.SOUTHERN))
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
if(params.diff_corr.isON)
    orbit = differential_correction(orbit.y0(1:6), cr3bp, orbit, params, cst);
else
    % Period (to be add into the .dat and directly interpolated):
    % Integration over one 1/2 orbit
    options = odeset('Events',@odezero_y,'Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
    [~,~,te] = ode45(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),[0 10],orbit.y0(1:6),options);
    orbit.T12 = te; %1/2 period
    orbit.T = 2*te; %period
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
if(params.plot.halo_orbit == cst.TRUE) %plotting  
  orbit_plot(yv, orbit, params, cst);
end
end