%-------------------------------------------------------------------------%
% Initializes a halo orbit by its vertical extension Az
%-------------------------------------------------------------------------%
% @return the structure orbit
function orbit = init_vert_lyap_orbit(cr3bp, li, family, Azdim, cst)

%CR3BP
orbit.cr3bp = cr3bp;

%Family
orbit.family = family;

%Class
switch(family)
    case cst.orbit.family.NORTHERN
        orbit.m = 1;
    case cst.orbit.family.SOUTHERN
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
orbit.Azdim = Azdim;
orbit.Az = Azdim/cr3bp.L;
%Az estimate
orbit.Azdim_estimate = Azdim;
orbit.Az_estimate = Azdim/cr3bp.L;

%Li
orbit.li = li;

%Status
orbit.status = cst.orbit.EMPTY;

end