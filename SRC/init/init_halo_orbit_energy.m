%-------------------------------------------------------------------------%
% Initializes a halo orbit by its jacobi constant
%-------------------------------------------------------------------------%
% @return the structure orbit
function orbit = init_halo_orbit_energy(cr3bp, li, family, Cjac, cst)

%CR3BP
orbit.cr3bp = cr3bp;

%family
orbit.family = family;

%Class
switch(family)
    case cst.orbit.NORTHERN
        orbit.m = 1;
    case cst.orbit.SOUTHERN
        orbit.m = 3;
end

switch li.number
    case 1
        orbit.dm = 2-orbit.m;
    case 2
        orbit.dm = orbit.m-2;  %BEWARE: Northern and Southern classes are opposite to Class I and Class II formulation (see Richardson 1980)
    case 3
        orbit.dm = 2-orbit.m;
    otherwise
        orbit.dm = 0;  %NO dm if Li.number!=1,2,3
end

%Az
orbit.C = Cjac;
orbit.E = -0.5*Cjac;

%Li
orbit.li = li;

%Status
orbit.status = cst.orbit.EMPTY;

end