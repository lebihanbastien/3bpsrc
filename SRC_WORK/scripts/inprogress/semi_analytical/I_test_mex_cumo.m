%--------------------------------------------------------------------------
% Example n°1: this matlab file makes use of the richardson third order
% approximation to build:
% - an EML2 planar lyapunov orbit (~8000 km of radius)
% - an EML2 halo orbit (~10 000km of vertical extension)
% - an EML2 vertical lyapunov orbit (~30 000km of vertical extension)
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Orbit init & computation for a planar lyapunov orbit
%Initialization
plyap = init_orbit(cr3bp, ...                     % Parent CR3BP
    cr3bp.l1, ...                  % Parent libration point
    cst.orbit.type.PLYAP, ...      % Planar lyapunov orbit
    cst.orbit.family.PLANAR, ...   % Planar class (useless here, since it is a planar lyapunov orbit
    8000, ...                      % Of Ax extension ~ 8000 km
    cst);                          % Numerical constants

%Computation
plyap = orbit_computation(cr3bp, plyap, default, cst);

%% Same for a halo orbit
%Initialization
halo = init_orbit(cr3bp, ...                     % Parent CR3BP
    cr3bp.l1, ...                  % Parent libration point
    cst.orbit.type.HALO, ...       % HALO orbit
    cst.orbit.family.NORTHERN, ... % Northern class
    1000, ...                     % Of vertical extension ~ 10000 km
    cst);                          % Numerical constants

%Computation
halo = orbit_computation(cr3bp, halo, default, cst);


%% Test of cmo.cpp

% For CRTBP EML1: st0_max = 0.25 (order 40)
% For CRTBP EML2: st0_max = 0.17 (order 40)

% @TODO: 1. CRTBP compute up to order 40 again (currently 30)
%        2. Check the order for QBCP (currently unknown, probably 20).
%        3. Get the manifold ouput as a 3-dim matrix? (gives the manifold leg output on
%        each point given in tspan (variable size) tspan = [t0 t1 t2...
%        tf]. Look for mxCreateNumericArray(dims[])
%        4. Change the computation: if output is just tf, yf, do NOT get
%        the results on a grid! Extra work for nothing!
%        5. Do the SEML case!

%----------
% Parameters
%----------
% Initial conditions
st0_max = 2.0;
st0     = st0_max*ones(1,4);
st0(5)  = 0.0;

% Libration point
li = cr3bp.l2;
% Order of the expansions
order = 12;
% Time span
tspan = [0 10*pi];
% Time on the manifold
tman = 5*pi;

% Type of coordinate system as in the output
outputType = cst.coord.VSYS;
% Type of model
model = 1; %CRTBP
% Type of framework
fwrk = 0; %Earth-Moon
% Man type
manType = 2;

%----------
%Get initial conditions from mex file
%----------
for k = 204:2:204
    %----------
    %Computation
    %----------
    [to, yo, tm, ym] = cumo(st0, tspan, tman, order, li.number, outputType, model, fwrk, manType, k);
    %----------
    %Plot
    %----------
    Lf = cr3bp.L*10^(-3);
    index = 3;
    figure(4);
    hold on
    axis equal
    grid on
    %Orbit branch
    plot3(yo(:,1)*Lf, yo(:,2)*Lf, yo(:,3)*Lf, 'Color', rgb('dark green'), 'LineWidth', 2);
    %Manifold branch
    plot3(ym(:,1)*Lf, ym(:,2)*Lf, ym(:,3)*Lf, 'Color', rgb('dark red'), 'LineWidth', 2);
end

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
    figure(4);
    view([-47 28]);
end