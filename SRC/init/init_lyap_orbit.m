%-------------------------------------------------------------------------%
% Initializes a planar lyapunov orbit by its x extension Ax
%-------------------------------------------------------------------------%
% @param cr3bp the structure containing the parent CR3BP
% @param li the number of the parent libration point (1 or 2)
% @param Axdim the desired approximated x-extension of the orbit (in km)
% @param cst the structure containing the numerical constants
% @return the structure orbit
function orbit = init_lyap_orbit(cr3bp, li, Axdim, cst)

%CR3BP
orbit.cr3bp = cr3bp;

%Ax
orbit.Axdim = Axdim;
orbit.Ax = Axdim/cr3bp.L;
%Ax estimate
orbit.Axdim_estimate = Axdim;
orbit.Ax_estimate = Axdim/cr3bp.L;

%Li
orbit.li = li;

%Status
orbit.status = cst.orbit.EMPTY;

end