function  ysg  = synodical2selenographic(ysyn, cr3bp)
%
% YSG  = SYNODICAL2SELENOGRAPHIC(YSYN, CR3BP) computes the coordinates YSYN
% into selenographic coordinates.
%
% BLB 2016

% Z-rotation with respect to synodical frame 
R01  = [-1 0 0 ; 0 -1 0 ; 0 0 1];
% Y-rotation with respect to synodical frame
tilt = deg2rad(6.687);
R12  = [cos(tilt) 0 -sin(tilt) ; 0 1 0 ; sin(tilt) 0  cos(tilt)];
R02  = R01*R12;

% Synodical state centered at the Moon
yt = ysyn - [cr3bp.m2.pos 0 0 0]';
% Result
ysg = zeros(6,1);

%Position
ysg(1:3) = R02\yt(1:3);
%Velocity
ysg(4:6) = R02\yt(4:6);

end

