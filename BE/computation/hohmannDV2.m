function [ dv2 ] = hohmannDV2( r1, r2, mu2BP )
% DV to leave the elliptical orbit of perigee r1 and apogee r2 and to enter
% the circular orbit of radius r2;
dv2 = sqrt(mu2BP/r2)*abs(1 - sqrt(2*r1/(r1+r2)));
end

