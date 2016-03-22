%--------------------------------------------------------------------------
% Continuation example n°3: Continuation procedure to produce a discrete set within 
% the family of northern orbits of EML2
% WARNING: very sensitive, works only for small Az amplitudes (contrary to
% planar lyapunov case). Complete continuations are more complex than that
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Change from default parameters
default.plot.orbit = false;
default.plot.XY = false;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Orbit
orbit_1 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.VLYAP, cst.orbit.family.NORTHERN, 16000, cst);
orbit_2 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.VLYAP, cst.orbit.family.NORTHERN, 20000, cst);

%% Orbit computation
orbit_1 = orbit_computation(cr3bp, orbit_1, default, cst);
orbit_2 = orbit_computation(cr3bp, orbit_2, default, cst);

%% Continuation
%Waitbar
h = waitbar(0,'Computation in progress...');

%Loop
maxIter = 500;
for i = 1:maxIter
    disp(i);
    yv = orbit_2.y0 + (orbit_2.y0 - orbit_1.y0);
    orbit_1 = orbit_2;
    
    if(mod(i,20) ==0)
        default.plot.orbit = true;
        orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yv, cst);
        default.plot.orbit = false;
    else
        orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yv, cst);
    end
    
    waitbar(i / maxIter);
    
end
close(h)

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end
