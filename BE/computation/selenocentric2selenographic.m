function ysg  = selenocentric2selenographic(ysc, t)
%
% YSG  = SELENOCENTRIC2SELENOGRAPHIC(YSC, CR3BP) computes the coordinates 
% YSC into selenographic coordinates.
%
% BLB 2016

%Rotation matrix
R = [cos(t) sin(t) 0 ; -sin(t) cos(t) 0 ; 0 0 1];
Rdot = [-sin(t) cos(t) 0 ; -cos(t) -sin(t) 0 ; 0 0 0];

ysg = zeros(6,1);
% Position
ysg(1:3) = R*ysc(1:3);

% Velocity
ysg(4:6) = R*ysc(4:6) + Rdot*ysc(1:3);

end
