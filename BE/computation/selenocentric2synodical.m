function  ysyn  = selenocentric2synodical(ysc, cr3bp, t)
%
% YSYN  = SELENOCENTRIC2SYNODICAL(YSC) computes the coordinates YSYN
% into synodical coordinates, from selenocentric coordinates.
%
% BLB 2016

%--------------------------------------------------------------------------
% Matrices
%--------------------------------------------------------------------------
% Z-rotation between the inertial frame and the synodical frame.
R = [cos(t) sin(t) 0 ; -sin(t) cos(t) 0 ; 0 0 1];
Rdot = [-sin(t) cos(t) 0 ; -cos(t) -sin(t) 0 ; 0 0 0];

% Z-rotation with respect to inertial frame to get opposite x/y
R01  = [-1 0 0 ; 0 -1 0 ; 0 0 1];

% Y-rotation with respect to inertial frame, to get a 6.687° tilt
tilt = deg2rad(6.687);
R12  = [cos(tilt) 0 -sin(tilt) ; 0 1 0 ; sin(tilt) 0  cos(tilt)];
R02  = R01*R12;

%--------------------------------------------------------------------------
% Inertial state centered at the Moon
%--------------------------------------------------------------------------
yin = zeros(6,1);

%Position
yin(1:3) = R02*ysc(1:3);
%Velocity
yin(4:6) = R02*ysc(4:6);

%--------------------------------------------------------------------------
% Result: synodical state
%--------------------------------------------------------------------------
ysynmc = zeros(6,1);

% Position
ysynmc(1:3) = R*yin(1:3);

% Velocity
ysynmc(4:6) = R*yin(4:6) + Rdot*yin(1:3);

% Centered at the Moon
ysyn = ysynmc + [cr3bp.m2.pos 0 0 0]';

end