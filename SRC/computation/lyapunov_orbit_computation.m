%-------------------------------------------------------------------------%
% Lyapunov orbits:
% Computation of orbit.y0 and various orbital parameters from a third-order 
% guess 
% WARNING: the use of third-order approximation should be limited to
% Ax < 20000km
%-------------------------------------------------------------------------%
% @param cr3bp the structure containing the parent CR3BP
% @param orbit the structure containing the parent orbit
% @param params the structure containing the comptuation parameters
% @param cst the structure containing the numerical constants
% @return the udpated structure orbit
function orbit = lyapunov_orbit_computation(cr3bp, orbit, params, cst)
%Third order approximation
yv_to =  third_order_lyap_orbit(orbit, 0.0);

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
    
%     %Second primary
%     Rm2 = cr3bp.m2.Req/cr3bp.L;
%     VTheta = 0:0.01:2*pi;
%     X_L = 1-cr3bp.mu + Rm2 *cos(VTheta);
%     Y_L = 0 + Rm2*sin(VTheta);
%     Z_L = 0 + Rm2*cos(VTheta);
%     
%     %Libration point
%     Li = orbit.li.position;
%     
%     figure(1)
%     hold on
%     %Orbit
%     plot(yv(:,1)*cr3bp.L*10^(-3),yv(:,2)*cr3bp.L*10^(-3), 'b');
%     %Second primary
%     fill(X_L*cr3bp.L*10^(-3), Y_L*cr3bp.L*10^(-3),'k');
%     plot(Li(1)*cr3bp.L*10^(-3), Li(2)*cr3bp.L*10^(-3),'rx');
%     %Libration point
%     xlabel('X (x 10^3 km)')
%     ylabel('Y (x 10^3 km)')
%     axis equal
%     grid on
%     title('Planar Lyapunov in X-Y plane');
%     
%     if(~params.plot.XY_only)
%     
%     figure(2)
%     hold on
%     %Orbit
%     plot(yv(:,1)*cr3bp.L*10^(-3),yv(:,3)*cr3bp.L*10^(-3), 'b', 'LineWidth',1.5);
%     %Second primary
%     fill(X_L*cr3bp.L*10^(-3), Y_L*cr3bp.L*10^(-3),'k');
%     %Libration point
%     plot(Li(1)*cr3bp.L*10^(-3), Li(3)*cr3bp.L*10^(-3),'rx')
%     xlabel('X (x 10^3 km)')
%     ylabel('Z (x 10^3 km)')
%     axis equal
%     grid on
%     title('Planar Lyapunov  in X-Z plane');
%     
%     
%     figure(3)
%     hold on
%     %Orbit
%     plot(yv(:,2)*cr3bp.L*10^(-3),yv(:,3)*cr3bp.L*10^(-3), 'b', 'LineWidth',1.5);
%     %Second primary
%     fill(Y_L*cr3bp.L*10^(-3), Z_L*cr3bp.L*10^(-3),'k');
%     xlabel('Y (x 10^3 km)')
%     ylabel('Z (x 10^3 km)')
%     axis equal
%     grid on
%     title('Planar Lyapunov orbit in Y-Z plane');
%     
%     end
%     
%     if(params.plot.TD)  %3D plot
%         
%         %Second primary
%         [X_E3D, Y_E3D, Z_E3D] = sphere;
%         X_E3D = 1-cr3bp.mu + Rm2*X_E3D;
%         Y_E3D = 0          + Rm2*Y_E3D;
%         Z_E3D = 0          + Rm2*Z_E3D;
%         
%         figure(4)
%         hold on
%         axis equal
%         grid on
%         xlabel('X (x 10^3 km)')
%         ylabel('Y (x 10^3 km)')
%         zlabel('Z (x 10^3 km)')
%         
%         %Second primary
%         surf(X_E3D*cr3bp.L*10^(-3), Y_E3D*cr3bp.L*10^(-3), Z_E3D*cr3bp.L*10^(-3), 'FaceColor', [23 153 179]./255, 'FaceLighting', 'none', 'EdgeColor', [9 63 87]./255);
%         %Orbit
%         plot3(yv(:,1)*cr3bp.L*10^(-3),yv(:,2)*cr3bp.L*10^(-3),yv(:,3)*cr3bp.L*10^(-3), 'b', 'LineWidth', 1.5);
%         %Libration point
%         plot3(Li(1)*cr3bp.L*10^(-3), Li(2)*cr3bp.L*10^(-3), Li(3)*cr3bp.L*10^(-3), 'rx', 'MarkerSize', 5);
%         %Set the camera
%         view([37.5 30])
%         %Zlimit
%         zlim([-50 50]);
%     end
%     
% end
end