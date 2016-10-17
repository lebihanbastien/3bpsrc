%--------------------------------------------------------------------------
% Example: this matlab file build some examples of orbits around a given
% libration point and save the corresponding figure. A manifold of one of
% the orbits is computed. note that the plots of the orbits are regrouped
% at the end.
%
% Author: BLB
% Year: 2016
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);
% Libration point
li = cr3bp.l1;

%% Plot init (figure 4)
initplot3D(4, cr3bp, default, li)
Lf = cr3bp.L;
figure(4);
hold on
axis equal
grid on

%% Change some default parameters
default.plot.orbit = false;
linewidth = 1.2;

%% Halo orbit
%--------------------------------------------------------------------------
%Initialization & computation
%--------------------------------------------------------------------------
halo = init_orbit(cr3bp,  li, cst.orbit.type.HALO,  cst.orbit.family.NORTHERN, 15000, cst);
halo = orbit_computation(cr3bp, halo, default, cst);

%% Planar Lyapunov orbit
%--------------------------------------------------------------------------
%Initialization & computation
%--------------------------------------------------------------------------
plyap = init_orbit(cr3bp, li, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, 8000, cst);
plyap = orbit_computation(cr3bp, plyap, default, cst);



%% Vertical Lyapunov orbit

%--------------------------------------------------------------------------
%Initialization & computation
%--------------------------------------------------------------------------
vlyap = init_orbit(cr3bp, li, cst.orbit.type.VLYAP, cst.orbit.family.NORTHERN, 30000, cst); 
vlyap = orbit_computation(cr3bp, vlyap, default, cst);

%% Lissajous orbit at EML1

%--------------------------------------------------------------------------
% For EML1: st0_max = 0.25 (order 40)
% For EML2: st0_max = 0.17 (order 40) 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Additionnal parameters for the cmo routine
%--------------------------------------------------------------------------
% Order of the expansions
order = 20;
% Time span (integration)
tspan = [0 30*pi];
% Type of coordinate system as in the output
outputType = cst.coord.VSYS; 
% Type of model
model = cst.model.RTBP; %CRTBP
% Type of framework
fwrk = cst.fwrk.EM; %Earth-Moon

%--------------------------------------------------------------------------
% Initial conditions
%--------------------------------------------------------------------------
st0_max = 0.15;
st0 = st0_max*ones(1,4);

%--------------------------------------------------------------------------
%Computation with cmo
%--------------------------------------------------------------------------
[~, yo] = cmo(st0, tspan, order, li.number, outputType, model, fwrk);

%% Manifolds & plots

%--------------------------------------------------------------------------
% Choose the orbit
%--------------------------------------------------------------------------
orbit = plyap;

%--------------------------------------------------------------------------
% Manifold initialization
%--------------------------------------------------------------------------

% We define a stable manifold
manifold_branch_stable    = init_manifold_branch(cst.manifold.STABLE, ...
                                                 cst.manifold.EXTERIOR);
% We define an unstable manifold                                            
manifold_branch_unstable  = init_manifold_branch(cst.manifold.UNSTABLE,...
                                                 cst.manifold.EXTERIOR);
                                             
% We define an event structure that can trigger the end of the integration 
% along the manifold. Here, the integration will stop when the line that
% links the current state and the center of the Earth reaches a given angle
% with respect to the Earth-Moon line.
moon.event = init_event(cst.manifold.event.type.ANGLE_SECTION,...         %the event is triggered at a given angle...
                        deg2rad(120),...                                    %of a given value...
                        cst.manifold.event.isterminal.YES,...              %the trajectory stops at the first ocurrence...
                        cst.manifold.event.direction.ALL,...               %all direction are considered...
                        cr3bp.m2.pos, cst);                                %the center for the computation of the angle is the Earth

% Note: the event can be directly incorporated to the manifold during its
% initialization (see init_manifold_branch.m).

%--------------------------------------------------------------------------
% Manifold computation                                             
%--------------------------------------------------------------------------
% Integration duration
t = 5;

%%Stable
% for theta = 0:0.05:1 %position on the orbit
%     manifold_branch_stable = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable, theta, t, default, cst, moon.event);
% end

%Unstable, with an event triggered.
for theta = 0:0.01:1 %position on the orbit
    manifold_branch_unstable = manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable, theta, t, default, cst, moon.event);
end

%% Plots of the orbits

plot3(halo.yv(:,1)*Lf, halo.yv(:,2)*Lf, halo.yv(:,3)*Lf    , 'LineWidth', linewidth, 'Color', rgb('gray'));
plot3(vlyap.yv(:,1)*Lf, vlyap.yv(:,2)*Lf, vlyap.yv(:,3)*Lf , 'LineWidth', linewidth, 'Color', rgb('gray'));
plot3(plyap.yv(:,1)*Lf, plyap.yv(:,2)*Lf, plyap.yv(:,3)*Lf ,'LineWidth', 1.5);%, 'LineWidth', linewidth, 'Color', rgb('gray'));
plot3(yo(:,1)*Lf, yo(:,2)*Lf, yo(:,3)*Lf                   , 'LineWidth', linewidth, 'Color', rgb('gray'));

%% Print
figtoprint(gcf, 'plot/earth_moon_orbits_plyap_manifold', [-150 14]);