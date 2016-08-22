function  ysc  = selenographic2selenocentric(ysg, t)
%
% YSC  = SELENOGRAPHIC2SELENOCENTRIC(YSG, CR3BP) computes the coordinates 
% YSG into selenocentric coordinates.
%
% BLB 2016

%Rotation matrix
R = [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1];
Rdot = [-sin(t) -cos(t) 0 ; cos(t) -sin(t) 0 ; 0 0 0];

% Position
ysc(1:3) = R*ysg(1:3);

% Velocity
ysc(4:6) = R*ysg(4:6) + Rdot*ysg(1:3);

end
