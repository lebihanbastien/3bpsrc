%--------------------------------------------------------------------------
% Example: this matlab file build some examples of orbits around a given
% libration point and save the corresponding figure. A manifold of the
% lissajous orbit is computed. Careful: do not change much until cumovec
% is updated (september 2016).
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

%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
plot3(halo.yv(:,1)*Lf, halo.yv(:,2)*Lf, halo.yv(:,3)*Lf, 'LineWidth', linewidth, 'Color', rgb('gray'));

%% Planar Lyapunov orbit
%--------------------------------------------------------------------------
%Initialization & computation
%--------------------------------------------------------------------------
plyap = init_orbit(cr3bp, li, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, 8000, cst);
plyap = orbit_computation(cr3bp, plyap, default, cst);

%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
plot3(plyap.yv(:,1)*Lf, plyap.yv(:,2)*Lf, plyap.yv(:,3)*Lf, 'LineWidth', linewidth, 'Color', rgb('gray'));

%% Vertical Lyapunov orbit

%--------------------------------------------------------------------------
%Initialization & computation
%--------------------------------------------------------------------------
vlyap = init_orbit(cr3bp, li, cst.orbit.type.VLYAP, cst.orbit.family.NORTHERN, 30000, cst); 
vlyap = orbit_computation(cr3bp, vlyap, default, cst);

%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
plot3(vlyap.yv(:,1)*Lf, vlyap.yv(:,2)*Lf, vlyap.yv(:,3)*Lf, 'LineWidth', linewidth, 'Color', rgb('gray'));

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
outputType = cst.coord2.VSYS; 
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


%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
plot3(yo(:,1)*Lf, yo(:,2)*Lf, yo(:,3)*Lf, 'LineWidth', linewidth);%, 'Color', rgb('gray'));

%--------------------------------------------------------------------------
% Manifold
%--------------------------------------------------------------------------
% Manifold type
manType    = cst.mantype.MAN_CENTER_U;
st0(5)     = 0.0;
% Time on the manifold (in Sun period)
tman = 1.0;

%Steps along the orbit
kmin  = 0;
kmax  = 150;
kstep = 2;

%--------------------------------------------------------------------------
%Computation
%--------------------------------------------------------------------------
[to, yo, tm, ym] = cumovec(st0, tspan, tman, order, li.number, outputType, model, fwrk, manType, kmin, kmax, kstep);

%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
for ii = 1:size(ym,3)
    plot3(ym(:,1,ii)*Lf, ym(:,2,ii)*Lf, ym(:,3,ii)*Lf, 'LineWidth', linewidth, 'Color', rgb('light red'));
end

%% Print
figtoprint(gcf, 'plot/earth_moon_orbits_lissajous_manifold', [-150 14]);