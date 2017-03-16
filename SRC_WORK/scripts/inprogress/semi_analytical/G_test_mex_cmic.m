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
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.PLYAP, ...      % Planar lyapunov orbit
    cst.orbit.family.PLANAR, ...   % Planar class (useless here, since it is a planar lyapunov orbit
    8000, ...                      % Of Ax extension ~ 8000 km
    cst);                          % Numerical constants

%Computation
plyap = orbit_computation(cr3bp, plyap, default, cst);

%% Same for a halo orbit
%Initialization
halo = init_orbit(cr3bp, ...                     % Parent CR3BP
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.HALO, ...       % HALO orbit
    cst.orbit.family.NORTHERN, ... % Northern class
    15000, ...                     % Of vertical extension ~ 10000 km
    cst);                          % Numerical constants

%Computation
halo = orbit_computation(cr3bp, halo, default, cst);


%% Test of cmic.cpp
%--------------------------------------------------------------------------
% @TODO:
%   1. Mettre en place un integrateur par projection sous C++ utilisant
%   directement RCMtoNC (pas de by TFC encore)
%   2. Le faire à pas de projection variable. Si l'erreur de projection
%   depasse un certain niveau, on relance la simulation au pas precedent
%   avec un delai de projection plus court. Si le delai de projection passe
%   sous une valeur fournie par l'utilisateur on arrete le calcul avec une
%   erreur.
%--------------------------------------------------------------------------
st0 = [0.0 0.3 0.0 0.3];
index = 1;
outputType = cst.coord2.VSYS;
li = cr3bp.l2;
t0 = 0.0;

for order = 10:10:40
    
    %----------
    %Get initial conditions from mex file
    %----------
    y0  = cmic(st0, t0, order, li.number, outputType);
    
    %----------
    %Integration
    %----------
    switch(outputType)
        case cst.coord2.NC
            options = odeset('Reltol', default.ode45.RelTol, 'Abstol', default.ode45.AbsTol);
            [~,yv] = ode45(@(t,y)cr3bp_mn_6(t,y, cr3bp.mu, li.c1, li.gamma_i) ,[0 5], y0, options);
        case cst.coord2.VSYS
            [~, yf, ~, yv] = ode78_cr3bp([0 5], y0, cr3bp.mu);
    end
    
    
    %----------
    %Factor
    %----------
    Lf = cr3bp.L*10^(-3);
    
    %----------
    %Plot
    %----------
    figure(4);
    hold on
    axis equal
    grid on
    %Manifold branch
    p(index) = plot3(yv(:,1)*Lf, yv(:,2)*Lf, yv(:,3)*Lf, 'Color', rgb(index), 'LineWidth', 2);
    legend(p,'10', '20', '30', '40');
    index=index+1;
end

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
    figure(4);
    view([-47 28]);
end