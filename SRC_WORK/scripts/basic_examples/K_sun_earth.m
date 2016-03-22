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

%% Inner changes from default parameters

% 1. Decomment the next line to use only MATLAB routines (very slow!)
%--------------------------------------------------------------------------
%default.computation.type = cst.computation.MATLAB;

% 2. Decomment the next lines to change the absolute and relative precision during integration with MATLAB routines  (ode45)
%--------------------------------------------------------------------------
default.ode45.AbsTol = 1e-12;
default.ode45.RelTol = 1e-12;

% 3. Decomment the next lines to change the plotting options
%--------------------------------------------------------------------------
default.plot.XY        = false; %plot also the results in X-Z plane

% 4. See parameters_default_init.m to see other options
%--------------------------------------------------------------------------

%% Environment init
cr3bp = init_CR3BP('SUN', 'EARTH', default);


%% Same for a halo orbit
%Initialization
isOrbitOnly = 0;
for kalt = 120000:10000:170000
halo = init_orbit(cr3bp, ...       % Parent CR3BP
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.HALO, ...       % HALO orbit
    cst.orbit.family.NORTHERN, ... % Northern class
    kalt, ...                      % Of vertical extension ~ 10000 km
    cst);                          % Numerical constants

%Computation
halo = orbit_computation(cr3bp, halo, default, cst, isOrbitOnly);
if(~isOrbitOnly) 
    isOrbitOnly = 1;
end
end
%Initialization
isOrbitOnly = 0;
for kalt = 120000:10000:170000
    halo = init_orbit(cr3bp, ...                     % Parent CR3BP
        cr3bp.l1, ...                  % Parent libration point
        cst.orbit.type.HALO, ...       % HALO orbit
        cst.orbit.family.NORTHERN, ... % Northern class
        kalt, ...                         % Of vertical extension ~ 10000 km
        cst);                          % Numerical constants
    
    %Computation
    halo = orbit_computation(cr3bp, halo, default, cst, isOrbitOnly);
    if(~isOrbitOnly) 
    isOrbitOnly = 1;
end
end
%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
    figure(4);
    view([47 28]);
end

%% Go on with matlab 2015a
h = figure(4);
set(h, 'Color',[1 1 1]);
set(h, 'Position', [100, 100, 1050, 900]);
zoom(2)
hold on;
grid off;
axis off;

%% Export
% h = prepareFig(h);
% exportFig(h,'CRTBP');