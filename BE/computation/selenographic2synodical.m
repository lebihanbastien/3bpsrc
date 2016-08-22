function  ysyn  = selenographic2synodical(ysg, cr3bp)
%
% YSYNMC  = SYNODICAL2SELENOGRAPHIC(YSG) computes the coordinates YSYN
% into synodical coordinates.
%
% BLB 2016

% Z-rotation with respect to synodical frame 
R01  = [-1 0 0 ; 0 -1 0 ; 0 0 1];
% Y-rotation with respect to synodical frame
tilt = deg2rad(6.687);
R12  = [cos(tilt) 0 -sin(tilt) ; 0 1 0 ; sin(tilt) 0  cos(tilt)];
R02  = R01*R12;

% Result
ysynmc = zeros(6,1);

%Position
ysynmc(1:3) = R02*ysg(1:3);
%Velocity
ysynmc(4:6) = R02*ysg(4:6);

%Centered at the Moon
ysyn = ysynmc + [cr3bp.m2.pos 0 0 0]';

end