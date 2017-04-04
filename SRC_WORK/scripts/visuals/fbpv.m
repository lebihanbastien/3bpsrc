%--------------------------------------------------------------------------
% Four-body problem introduction
%
% BLB 2016
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters

%--------------------------------------------------------------------------
% Decomment the next lines to change the plotting options
%--------------------------------------------------------------------------
default.plot.XY             = true; % plot also the results in X-Z plane
default.plot.XZ             = false; % plot also the results in X-Z plane
default.plot.YZ             = false; % plot also the results in Y-Z plane
default.plot.TD             = false; % plot also the results in Y-Z plane
default.plot.diff_corr      = false; % plot the differential correction steps
default.plot.firstPrimDisp  = true;  % plot the Earth

%--------------------------------------------------------------------------
% See init_parameters_default.m to see other options
%--------------------------------------------------------------------------

%% "Good" User parameters (provide a good connection) 
user.Az    = 10000;                 %size of the orbit (km)
user.tau   = 0.17;                  %position on the orbit (adim, between 0 and 1)
user.t0    = 8;                     %integration time on the manifold in cr3bp (adim)
user.tl    = 20;                    %integration time on the rest of the trajectory (adim)
user.initSunPos = degtorad(100);    %initial position of the Sun (rad)

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);
sem   = init_CR3BP('SUN', 'EARTH_AND_MOON', default);

%% Plot init

default.plot.firstPrimDisp  = false; %don't plot the Sun
initplot2D(5, sem, default, sem.l1, 1, 2);
default.plot.firstPrimDisp  = true; %back to true

%% halo orbit
%--------------------------------------------------------------------------
%Initialization
%--------------------------------------------------------------------------
halo = init_orbit(cr3bp, cr3bp.l2, cst.orbit.type.HALO, cst.orbit.family.NORTHERN, user.Az , cst);

%--------------------------------------------------------------------------
%Computation
%--------------------------------------------------------------------------
halo = orbit_computation(cr3bp, halo, default, cst);

%% Manifold leg: the good one
%--------------------------------------------------------------------------
% We define an unstable manifold
%--------------------------------------------------------------------------
mbu  = init_manifold_branch(cst.manifold.UNSTABLE, cst.manifold.EXTERIOR);

%--------------------------------------------------------------------------
% Event at the earth
%--------------------------------------------------------------------------
earth.event = init_event(cst.manifold.event.type.ANGLE_SECTION,...         %the event is triggered at a given angle...
                        deg2rad(150),...                                   %of a given value...
                        cst.manifold.event.isterminal.YES,...              %the trajectory stops at the first ocurrence...
                        cst.manifold.event.direction.ALL,...               %all direction are considered...
                        cr3bp.m1.pos, cst);                                %the center for the computation of the angle is the Earth
                    
%--------------------------------------------------------------------------
% Good Parameters
%--------------------------------------------------------------------------
% Position on the orbit (between 0 and 1)
tau = user.tau;
%Computation
mbu = manifold_branch_computation(cr3bp, halo, mbu, tau, user.t0, default, cst, earth.event);
te = mbu.te;

%--------------------------------------------------------------------------
% Loop
%--------------------------------------------------------------------------
for tau = 0:0.05:0.95
    %Computation 
    mbu = manifold_branch_computation(cr3bp, halo, mbu, tau, user.t0, default, cst, earth.event);
    mbu.ttraj = mbu.ttraj + (te - mbu.te);
    
    %Plot in the SEM framework
    [~, yv_sem ] = bcp_em_to_sem(mbu.ttraj, mbu.ytraj,  user.initSunPos, cr3bp, sem, cst);
    
    figure(5);
    plot(yv_sem(:,1)*sem.L, yv_sem(:,2)*sem.L, 'Color', rgb('light red'), 'LineWidth', 1.2);
end

%--------------------------------------------------------------------------
% Good Parameters
%--------------------------------------------------------------------------
% Position on the orbit (between 0 and 1)
tau = user.tau;
%Computation
mbu = manifold_branch_computation(cr3bp, halo, mbu, tau, user.t0, default, cst, earth.event);
%Plot in the SEM framework
[~, yv_sem ] = bcp_em_to_sem(mbu.ttraj, mbu.ytraj,  user.initSunPos, cr3bp, sem, cst);
figure(5);
plot(yv_sem(:,1)*sem.L, yv_sem(:,2)*sem.L, 'Color', rgb('light red'), 'LineWidth', 1.2);

%% Replot the orbit
orbit_plot(halo, default);

%% In the 4bp

%--------------------------------------------------------------------------
% Interval of integration
%--------------------------------------------------------------------------
tspan = [mbu.te mbu.te+user.tl];

%--------------------------------------------------------------------------
% Moreover, events structure can be used to stop the integration at a given 
% time or to store some events along the trajectory
%--------------------------------------------------------------------------
% Example: a null flight path angle (velocity vector orthogonal to the
% Earth-Spacecraft line (at perigee or apogee).
earth.event = init_event(cst.manifold.event.type.FLIGHT_PATH_ANGLE,...     %the event is triggered when the flight path angle is...
                         0.0,...                                           %equal to zero...
                         cst.manifold.event.isterminal.YES,...             %the trajectory stops at the first ocurrence...
                         cst.manifold.event.direction.ALL,...              %all direction are considered...
                         cr3bp.m1.pos, cst);                               %the center for the computation of the angle is the Earth

%--------------------------------------------------------------------------
% Initial position of the Sun (in radians)
%--------------------------------------------------------------------------
initSunPos = user.initSunPos;   

%------------------------------------
% With ode78_bcp and event
%------------------------------------
earth.event.max_events = 2;
tic()
[~, yearc_bcp, tarc_bcp, yarc_bcp] = ode78_bcp_event(tspan, mbu.yve(1:6), cr3bp.mu, initSunPos, cst.sun.ms, cst.sun.as, cst.sun.omegaS, earth.event);
dtoc = toc();
fprintf('ode78_bcp integration in %5.5f s\n', dtoc);

% CAREFUL! Note that the outputs are inversed with respect to the MATLAB
% implementation. It was easier to do that way in compiled C.

%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------
% 2D
figure(1);
hold on
hp(2) = plot(yarc_bcp(:,1)*cr3bp.L, yarc_bcp(:,2)*cr3bp.L, 'Color', rgb('dark green'), 'LineWidth', 1.2);
if(exist('yearc_bcp', 'var'))
    plot(yearc_bcp(:,1)*cr3bp.L, yearc_bcp(:,2)*cr3bp.L, 'o', 'Color', rgb('dark green'), 'MarkerFaceColor', rgb('dark green'), 'MarkerSize', 3);
end

%In SEM framework
[~, yv_sem ] = bcp_em_to_sem(tarc_bcp, yarc_bcp,  user.initSunPos, cr3bp, sem, cst);
figure(5);
hold on
plot(yv_sem(:,1)*sem.L, yv_sem(:,2)*sem.L, 'Color', rgb('dark green'), 'LineWidth', 1.2);

%% Moon orbit, in SEM framework
%--------------------------------------------------------------------------
% Several orbits in a row
%--------------------------------------------------------------------------
ttot = [mbu.ttraj ; tarc_bcp];
yMoon = zeros(size(ttot, 1),6);

for p = 1:size(ttot, 1)
    yMoon(p,:) = [cr3bp.m2.pos 0 0 0];
end

[~, yMoon_sem] = bcp_em_to_sem(ttot, yMoon,  user.initSunPos, cr3bp, sem, cst);
figure(5);
hold on
plot(yMoon_sem(:,1)*sem.L, yMoon_sem(:,2)*sem.L, 'Color', rgb('black'), 'LineWidth', 1.2);


%% Print
figtoprint(figure(1), 'plot/coupled_crtbp_em', [0 90]);
figtoprint(figure(5), 'plot/coupled_crtbp_sem', [0 90]);

%% Halo orbit in SEM framework (not really necessary and does not work...)
% n = floor((ttot(end) - ttot(1))/halo.T);
% m = size(ttot, 1);
% yvc = halo.y0';
% tvc = ttot(1);
% for p = 1:n
%     ind = (p-1)*floor(m/n)+1; 
%     [~, ~, tv, yv] = ode78_cr3bp([ttot(ind) ttot(ind)+halo.T], halo.y0, cr3bp.mu);
%     [~, yv_sem] = bcp_em_to_sem(tv, yv(:,1:6),  user.initSunPos, cr3bp, sem, cst);
%     figure(5);
%     hold on
%     plot(yv_sem(:,1)*sem.L, yv_sem(:,2)*sem.L, 'Color', rgb('magenta'), 'LineWidth', 1.2);
% end


