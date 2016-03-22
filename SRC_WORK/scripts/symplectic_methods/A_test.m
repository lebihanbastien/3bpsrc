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
default.computation.type = cst.computation.MATLAB;  

% 2. Decomment the next lines to change the absolute and relative precision during integration with MATLAB routines  (ode45)
%--------------------------------------------------------------------------
default.ode45.AbsTol = 1e-10;
default.ode45.RelTol = 1e-10;

% 3. Decomment the next lines to change the plotting options
%--------------------------------------------------------------------------
%default.plot.XZ        = true; %plot also the results in X-Z plane
%default.plot.YZ        = true; %plot also the results in Y-Z plane
%default.plot.diff_corr = true; %plot the differential correction steps

% 4. See parameters_default_init.m to see other options
%--------------------------------------------------------------------------

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

%% Test of Euler B

% Initial conditions in canonical form
z0(1:3) = plyap.y0(1:3);             %position is unchanged
z0(4)   = plyap.y0(4) - plyap.y0(2); %px = xdot - y
z0(5)   = plyap.y0(5) + plyap.y0(1); %py = ydot + x
z0(6)   = plyap.y0(6);               %pz = zdot

% Euler-B scheme
dt = 0.1;
vS = floor(plyap.T/dt);

tvec = zeros(vS, 1);
zvec = zeros(vS, 6);

% Initial conditions
zvec(1,:) = z0;
tvec(1)   = 0;

%%

% Integration
for n = 1:vS-1
    
    tvec(n+1) = tvec(n)+dt; 
    
    %v^(n+1) = v^n - dt*Phi'(q^(n))
    %Only the vfn(2) is used
    vfn = LennardJonesVF(tvec(n), zvec(n,:));
    zvec(n+1,2) = zvec(n,2) + dt*vfn(2);
    
    
    % q^(n+1) = q^n + dt*p^(n+1)
    zvec(n+1,1) = zvec(n,1) + dt*zvec(n+1,2);
    
end
    
    
    
    
end

