%-------------------------------------------------------------------------%
% orbit_computation(cr3bp, orbit, params, cst)
%
% Computation of the initial conditions orbit.y0 and various orbital 
% parameters from a third-order guess.
%
% WARNING: the orbit structure should be properly initialized via the
% routine init_orbit.
%-------------------------------------------------------------------------%
% Inputs:
% 1. cr3bp  the structure containing the parent CR3BP
% 2. orbit  the structure containing the orbit to update
% 3. params the structure containing the computation parameters
% 4. cst the structure containing the numerical constants
%
% Outputs:
% 1. the updated orbit structure 
%
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function orbit = orbit_computation(cr3bp, orbit, params, cst, isOrbitOnly)
               
%-------------------------------------------------------------------------%
% Initialization from third order approximation
%-------------------------------------------------------------------------%                 
yv_guess =  third_order_orbit(orbit, 0.0, cst);

%-------------------------------------------------------------------------%
% Refinement
%-------------------------------------------------------------------------%     
orbit = orbit_refinement(cr3bp, orbit, params, yv_guess, cst, isOrbitOnly);

end