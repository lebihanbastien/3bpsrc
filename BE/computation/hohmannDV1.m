function [ dv1 ] = hohmannDV1( r1, r2, mu2BP )
% DV to enter the elliptical orbit of perigee r1 and apogee r2
dv1 = sqrt(mu2BP/r1)*abs(sqrt(2*r2/(r1+r2)) - 1);
end

