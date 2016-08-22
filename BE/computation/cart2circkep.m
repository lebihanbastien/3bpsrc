function [ a, e, I, omega, Omega ] = cart2circkep( ysc )
% Cartesian to Keplerian coordinates for a circular orbit only.
%
% BLB 2016

r = ysc(1:3);     %position in selenocentric frame
v = ysc(4:6);     %velocity in selenocentric frame
h =  cross(r, v); %Angular momentum vector

% Keplerian elements
Omega = atan2(h(1), -h(2));
I = atan2(norm(h(1:2)), h(3));
e = 0.0;
a = norm(r);
p = r2p(r', I, Omega);
omega = atan2(p(2), p(1));

end

