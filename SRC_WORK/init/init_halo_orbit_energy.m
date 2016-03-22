%-------------------------------------------------------------------------%
% init_halo_orbit_energy(cr3bp, li, family, Cjac, cst).
%
% Initializes a halo orbit via its Jacobi constant Cjac. Alternative to the
% routine init_orbit(cr3bp, li, type, Afdim, family, cst), which can
% initialize a halo orbit by an estimate of its vertical extension Az.
%-------------------------------------------------------------------------%
% Inputs:
% 
% 1. cr3bp the structure containing the parent CR3BP
% 2. li the structure containing the parent libration point.
% 3. family: the orbit family: either NORTHERN or SOUTHER for halo and
%    vertical lyapunov orbits. Always PLANAR for planar lyapunov orbits. 
%    Note that the family is forced to PLANAR for planar lyapunov orbits,
%    regardless of the input.
% 4. Azdim: The vertical extension of the orbit (in km)
% 5. cst the structure containing the numerical constants
%
% Outputs:
%
% 1. orbit the structure orbit, with the following fields:
%           - orbit.cr3bp       |
%           - orbit.li          |
%           - orbit.type        |
%           - orbit.C           |
%           - orbit.E           |
%           - orbit.family      |
%           - orbit.m           |
%           - orbit.dm          |
%           - orbit.status      |
%       
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function orbit = init_halo_orbit_energy(cr3bp, li, family, Cjac, cst)

%CR3BP
orbit.cr3bp = cr3bp;

%Li
orbit.li = li;

%Family
orbit.family = family;

%Type
orbit.type = cst.orbit.type.HALO;

%Class
switch(family)
    case cst.orbit.family.NORTHERN
        orbit.m = 1;
    case cst.orbit.family.SOUTHERN
        orbit.m = 3;
end

%dm parameter (see Richardson)
switch li.number
    case 1
        orbit.dm = 2-orbit.m;
    case 2
        orbit.dm = orbit.m-2;  %BEWARE: Northern and Southern classes are opposite to Class I and Class II formulation (see Richardson 1980)
    case 3
        orbit.dm = 2-orbit.m;
    otherwise
        orbit.dm = NaN;  %NO dm if Li.number!=1,2,3
end

%Energy: 
% - orbit.C = Jacobi constant
% - orbit.E = -1/2*orbit.C = energy of the orbit
orbit.C = Cjac;
orbit.E = -0.5*Cjac;

%Status
orbit.status = cst.orbit.EMPTY;

end