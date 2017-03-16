%--------------------------------------------------------------------------
% Four-body problem introduction
%
% BLB 2016
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters

% Decomment the next lines to change the plotting options
%--------------------------------------------------------------------------
default.plot.XZ             = false; % plot also the results in X-Z plane
default.plot.YZ             = false; % plot also the results in Y-Z plane
default.plot.diff_corr      = false; % plot the differential correction steps
default.plot.firstPrimDisp  = true;  % plot the Earth

% See parameters_default_init.m to see other options
%--------------------------------------------------------------------------

%% User parameter
user.Az    = 10000;                 %size of the orbit (km)
user.tau   = 0.01;                  %position on the orbit (adim, between 0 and 1)
user.t0    = 5;                     %integration time on the manifold in cr3bp (adim)
user.tl    = 10;                     %integration time on the rest of the trajectory (adim)
user.initSunPos = degtorad(100);    %initial position of the Sun (rad)

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% halo orbit
%Initialization
halo = init_orbit(cr3bp, ...       % Parent CR3BP
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.HALO, ...       % HALO orbit
    cst.orbit.family.NORTHERN, ... % Northern class
    user.Az , ...                  % Of vertical extension user.Az
    cst);                          % Numerical constants

%Computation
halo = orbit_computation(cr3bp, halo, default, cst);

%% Manifold leg
% We define an unstable manifold
manifold_branch_unstable  = init_manifold_branch(cst.manifold.UNSTABLE, cst.manifold.EXTERIOR);
% Integration duration
t0 = user.t0;

% Position on the orbit (between 0 and 1)
tau = user.tau;

manifold_branch_unstable = manifold_branch_computation(cr3bp, halo, manifold_branch_unstable, tau, t0, default, cst);

%% MATLAB ode integration options (important for precision!)
options = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol);

%% Rest of the trajectory can be computed in the CR3BP...

%--------------------------------------------------------------------------
% Interval of integration
%--------------------------------------------------------------------------
tspan = [user.t0 user.t0+user.tl];

%--------------------------------------------------------------------------
% Integration. 2 possibilities:
%       - With MATLAB ode113 routine
%       - With MEX file ode78_bcp (compiled C routine)
%--------------------------------------------------------------------------
% With ode113
%------------------------------------
tic()
[~, yarc_3bp] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),tspan, manifold_branch_unstable.yve, options);
dtoc = toc();
fprintf('ode113 integration in %5.5f s\n', dtoc);

%------------------------------------
% With ode78_cr3bp
%------------------------------------
% tic()
% [~, ~, ~, yarc_3bp] = ode78_cr3bp(tspan, manifold_branch_unstable.yv(1:6), cr3bp.mu);
% dtoc = toc();
% fprintf('ode78_cr3bp integration in %5.5f s\n', dtoc);

%--------------------------------------------------------------------------
% Moreover, events structure can be used to stop the integration at a given 
% time or to store some events along the trajectory
%--------------------------------------------------------------------------
% Example: a null flight path angle (velocity vector orthogonal to the
% Earth-Spacecraft line (at perigee or apogee).
earth.event = init_event(cst.manifold.event.type.FLIGHT_PATH_ANGLE,...     %the event is triggered when the flight path angle is...
                         0.0,...                                           %equal to zero...
                         cst.manifold.event.isterminal.NO,...             %the trajectory stops at the first ocurrence...
                         cst.manifold.event.direction.ALL,...              %all direction are considered...
                         cr3bp.m1.pos, cst);                               %the center for the computation of the angle is the Earth

%Include the event in the MATLAB computation is a bit tedious, it is easier
%to do with MEX routines (see below)
options_with_event = odeset('Events',@(t,y)odezero_flightpathangle(t,y,earth.event),...
                            'Reltol', default.ode113.RelTol,...
                            'Abstol', default.ode113.AbsTol);
                        
%------------------------------------
% Integration with ode113 and event
%------------------------------------
% tic()
% [~, yarc_3bp, ~, yearc_3bp] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),tspan, manifold_branch_unstable.yv, options_with_event);
% dtoc = toc();
% fprintf('ode113 integration in %5.5f s\n', dtoc);        

%--------------------------------------------------------------------------
% Again, MATLAB and MEX routines can be used.
% There is a small but useful difference between MATLAB and MEX
% implementation:
%  - With MATLAB (e.g. ode113, you can either get all events along the
%  trajectory OR stop at the very first one, setting event.isterminal=true.
% (this field is only used by MATLAB routines)
%  - With MEX files (e.g. ode78_bcp_event), you get a fixed number of
%  events along the trajectory, via the value  event.max_events (this field
%  is only used by MEX routines). If you want ALL events, just set
%  max_events to a very big values (the max allowed is 1000).
%--------------------------------------------------------------------------
earth.event.max_events = 2;
%------------------------------------
% With ode78_cr3bp_event and event
%------------------------------------
% tic()
% [~, yearc_3bp, ~, yarc_3bp] = ode78_cr3bp_event(tspan, manifold_branch_unstable.yv(1:6), cr3bp.mu, earth.event);
% dtoc = toc();
% fprintf('ode78_cr3bp_event integration in %5.5f s\n', dtoc);

% CAREFUL! Note that the outputs are inversed with respect to the MATLAB
% implementation. It was easier to do that way in compiled C.

%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------
%2D
figure(1);
hold on
hp(1) = plot(yarc_3bp(:,1)*cr3bp.L, yarc_3bp(:,2)*cr3bp.L, 'Color', rgb('dark blue'), 'LineWidth', 1.5);
if(exist('yearc_3bp', 'var'))
    plot(yearc_3bp(:,1)*cr3bp.L, yearc_3bp(:,2)*cr3bp.L, 'o', 'Color', rgb('dark blue'), 'MarkerFaceColor', rgb('dark green'), 'MarkerSize', 3);
end


%3D
figure(4);
hold on
ht(1) = plot3(yarc_3bp(:,1)*cr3bp.L, yarc_3bp(:,2)*cr3bp.L, yarc_3bp(:,3)*cr3bp.L, 'Color', rgb('dark blue'), 'LineWidth', 1.5);


%% ... Or in the 4bp

%--------------------------------------------------------------------------
% Initial position of the Sun (in radians)
%--------------------------------------------------------------------------
initSunPos = user.initSunPos;

%--------------------------------------------------------------------------
% Integration. Same remarks
%--------------------------------------------------------------------------
% With ode113
%------------------------------------
tic()
[~, yarc_bcp] = ode113(@(t,y)bcfbp_derivatives_6(t,y,cr3bp.mu, 0.2, cst.sun.ms, cst.sun.as, cst.sun.omegaS),tspan, manifold_branch_unstable.yve, options);
dtoc = toc();
fprintf('ode113 integration in %5.5f s\n', dtoc);

%------------------------------------
% With ode78_bcp
%------------------------------------
% tic()
% [~, ~, ~, yarc_bcp] = ode78_bcp(tspan, manifold_branch_unstable.yv(1:6), cr3bp.mu, initSunPos, cst.sun.ms, cst.sun.as, cst.sun.omegaS);
% dtoc = toc();
% fprintf('ode78_bcp integration in %5.5f s\n', dtoc);
 
%------------------------------------
% With ode113 and event
%------------------------------------
% tic()
[~, yarc_bcp, ~, yearc_bcp] = ode113(@(t,y)bcfbp_derivatives_6(t,y,cr3bp.mu, initSunPos, cst.sun.ms, cst.sun.as, cst.sun.omegaS),tspan, manifold_branch_unstable.yv, options_with_event);
% dtoc = toc();
% fprintf('ode113 integration in %5.5f s\n', dtoc);        

%------------------------------------
% With ode78_bcp and event
%------------------------------------
% earth.event.max_events = 3;
% tic()
[~, yearc_bcp, ~, yarc_bcp] = ode78_bcp_event(tspan, manifold_branch_unstable.yv(1:6), cr3bp.mu, initSunPos, cst.sun.ms, cst.sun.as, cst.sun.omegaS, earth.event);
% dtoc = toc();
% fprintf('ode78_bcp integration in %5.5f s\n', dtoc);

% CAREFUL! Note that the outputs are inversed with respect to the MATLAB
% implementation. It was easier to do that way in compiled C.

%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------
% 2D
figure(1);
hold on
hp(2) = plot(yarc_bcp(:,1)*cr3bp.L, yarc_bcp(:,2)*cr3bp.L, 'Color', rgb('dark green'), 'LineWidth', 1.5);
if(exist('yearc_bcp', 'var'))
    plot(yearc_bcp(:,1)*cr3bp.L, yearc_bcp(:,2)*cr3bp.L, 'o', 'Color', rgb('dark green'), 'MarkerFaceColor', rgb('dark green'), 'MarkerSize', 3);
end
legend(hp, {'CR3BP', 'BCP'});

% 3D
figure(4);
hold on
ht(2) = plot3(yarc_bcp(:,1)*cr3bp.L, yarc_bcp(:,2)*cr3bp.L, yarc_bcp(:,3)*cr3bp.L, 'Color', rgb('dark green'), 'LineWidth', 1.5);
legend(ht, {'CR3BP', 'BCP'});


%% Plot the 3bsoi sphere

%--------------------------------------------------------------------------
% 2D
%--------------------------------------------------------------------------
figure(1);
hold on
%Sphere
vt = 0:0.01:2*pi;
tbsoi      = zeros(2, size(vt,2));
tbsoi(1,:) = cst.env.em3bsoi * cos(vt);
tbsoi(2,:) = cst.env.em3bsoi * sin(vt);
%Position
vl = [1-cr3bp.mu  0  0]*cr3bp.L;
plot((vl(1) + tbsoi(1,:)),(vl(2) + tbsoi(2,:)), 'k--');

%--------------------------------------------------------------------------
% 3D
%--------------------------------------------------------------------------
figure(4);
hold on
[x, y, z] = sphere(16);
h = surf(cst.env.em3bsoi * x + vl(1), cst.env.em3bsoi * y, cst.env.em3bsoi * z);
set(h, 'FaceAlpha', 0.1, 'EdgeColor', 0.5*ones(3,1), 'FaceColor', 0.5*ones(3,1));
shading interp