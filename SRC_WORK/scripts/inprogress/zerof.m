function [ fz ] = zerof(x0, cr3bp, default, cst)
orbit = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, x0, cst);
orbit = orbit_computation(cr3bp, orbit, default, cst);
fz = orbit.T - 3.395596935953959;
fprintf('x0 = %5.15f \n', orbit.y0(1));
fprintf('vy0 = %5.15f \n', orbit.y0(5));
end

