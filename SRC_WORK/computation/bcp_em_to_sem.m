function [ tv_sem, yv_sem ] = bcp_em_to_sem( tv_em, yv_em,  initSunPos, cr3bp_em, cr3bp_sem, cst)
% Make the change of coordinates from the EM framework to the SEM framework
% in the Bicircular Problem.

%--------------------------------------------------------------------------
% Initialize some constants
%--------------------------------------------------------------------------
L_AB = 1/(cst.sun.as);   %distance ratio
T_AB = 1+cst.sun.omegaS; %time ratio

%--------------------------------------------------------------------------
% Init the outputs
%--------------------------------------------------------------------------
yv_sem  = zeros(size(yv_em, 1), size(yv_em, 2));
tv_sem  = zeros(size(tv_em, 1), size(tv_em, 2)); 
xIN_sem = zeros(6,1);
%--------------------------------------------------------------------------
% Loop on all the coordinates
%--------------------------------------------------------------------------
for p = 1:size(tv_em, 1)
    %----------------------------------------------------------------------
    % Computation
    %----------------------------------------------------------------------
    %Rotation matrix
    R_A = bcpRotMat(tv_em(p));
    
    %To inertial EM coordinates
    xIN_em = R_A*(yv_em(p,1:6)' - [cr3bp_em.m1.pos'; 0; 0; 0]);
    
    %To inertial SEM coordinates
    xIN_sem(1:3) = xIN_em(1:3)*L_AB;
    xIN_sem(4:6) = xIN_em(4:6)*L_AB/T_AB;
    
    %Current angle in SEM units
    theta_B = tv_em(p)*T_AB + initSunPos;
    
    %Rotation matrix, in SEM units
    R_B = bcpRotMat(theta_B);
    
    %----------------------------------------------------------------------
    % Outputs
    %----------------------------------------------------------------------
    %To SEM coordinates
    yv_sem(p,:) = R_B \ xIN_sem + [cr3bp_sem.m2.pos'; 0; 0; 0];
    
    % Time in SEM units
    tv_sem = tv_em(p)*T_AB;
end

%--------------------------------------------------------------------------
% Subroutines
%--------------------------------------------------------------------------
function RotMat = bcpRotMat(theta)
% Compute the rotation matrix associated to the angle theta
%
%   R =  | R11    0  |
%        | R21  R11  |
% with
%
%         | c -s 0 |          | -s -c 0 |
%   R11 = | s  c 0 |,   R21 = |  c -s 0 | 
%         | 0  0 1 |          |  0  0 0 |
% and
%       c = cos(theta), s = sin(theta)
%
% BLB 2016
RotMat11 = [+cos(theta) -sin(theta) 0; +sin(theta) +cos(theta) 0 ; 0 0 1];
RotMat21 = [-sin(theta) -cos(theta) 0; +cos(theta) -sin(theta) 0 ; 0 0 0];
RotMat   = [RotMat11 zeros(3) ; RotMat21 RotMat11];
end

end