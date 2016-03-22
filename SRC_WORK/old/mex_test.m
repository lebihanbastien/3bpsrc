%--------------------------------------------------------------------------
% Test of the mex routines
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%Number format
format long

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);



%% Orbit init & computation
lyap = init_lyap_orbit(cr3bp, cr3bp.l1, 8000 , cst);
lyap = lyapunov_orbit_computation(cr3bp, lyap, default, cst);

%% Computation over t1 = 1.0
options = odeset('Events',@odezero_y,'Reltol', default.ode45.RelTol, 'Abstol', default.ode45.AbsTol);

tic;
[te,ye] = ode45(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),[0 lyap.T12],lyap.y0(1:6),options);
T = toc;


disp('Via ode45, y(t1) = ')
disp(ye(end,:))
fprintf('computed in %f s\n', T)


%% Test of the mex function
tic;
[tf, yf] = ode78_cr3bp(0.0, lyap.T12, lyap.y0, 42, cr3bp.mu);
T = toc;

disp('Via mex function, y(t1) = ')
disp(yf(1:6))
fprintf('computed in %f s\n', T)


disp('Difference = ')
disp(ye(end,:)-yf(1:6))



%% Test of the mex function

%Event structure
val_par = init_event(cst.manifold.event.type.Y_SECTION,...
                     degtorad(0.0),...
                     cst.manifold.event.isterminal.YES,...
                     cst.manifold.event.direction.ALL,... 
                     cr3bp.m2.pos,...
                     cst);
    

tic;
[te, ye, tv, yv] = ode78_cr3bp_event(0.0, 10*lyap.T12, lyap.y0, 6, cr3bp.mu, val_par);
T = toc;



figure;
plot(yv(:,1), yv(:,2))
