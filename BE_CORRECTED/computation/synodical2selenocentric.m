function  ysc  = synodical2selenocentric(ysyn, cr3bp, t)
%
% YSC  = SYNODICAL2SELENOCENTRIC(YSYN, CR3BP) computes the coordinates YSYN
% into selenographic coordinates.
%
% BLB 2016

%--------------------------------------------------------------------------
% Matrices
%--------------------------------------------------------------------------
% Z-rotation between the inertial frame and the synodical frame.
R = [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1];
Rdot = [-sin(t) -cos(t) 0 ; cos(t) -sin(t) 0 ; 0 0 0];

% Z-rotation with respect to inertial frame to get opposite x/y
R01  = [-1 0 0 ; 0 -1 0 ; 0 0 1];

% Y-rotation with respect to inertial frame, to get a 6.687° tilt
tilt = deg2rad(6.687);
R12  = [cos(tilt) 0 -sin(tilt) ; 0 1 0 ; sin(tilt) 0  cos(tilt)];
R02  = R01*R12;

%--------------------------------------------------------------------------
% Synodical state centered at the Moon
%--------------------------------------------------------------------------
ysynmc = ysyn - [cr3bp.m2.pos 0 0 0]';

%--------------------------------------------------------------------------
% Inertial state centered at the Moon
%--------------------------------------------------------------------------
yin = zeros(6,1);
% Position
yin(1:3) = R*ysynmc(1:3);

% Velocity
yin(4:6) = R*ysynmc(4:6) + Rdot*ysynmc(1:3);

%--------------------------------------------------------------------------
% Result: selenocentric state
%--------------------------------------------------------------------------
%ysc = zeros(6,1);
%Position
ysc(1:3) = R02\yin(1:3);
%Velocity
ysc(4:6) = R02\yin(4:6);

end