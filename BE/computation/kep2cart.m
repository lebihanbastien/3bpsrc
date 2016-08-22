function [ yc ] = kep2cart(a, e, I, omega, Omega, M, mu2BP)

% Circular orbit: M  = E, for now
E = CalcEA(M,e,1e-10);

% Position in the orbital plane
q = a*[cos(E) - e ; sqrt(1-e*e)*sin(E) ; 0];

% Mean motion
n = sqrt(mu2BP/a^3);

% Velocity in the orbital plane;
qdot = n*a/(1 - e*cos(E))*[-sin(E) ; sqrt(1-e*e)*cos(E) ; 0];

% Back to inertial coordinates
yc = zeros(6,1);
yc(1:3) = R3(-Omega)*R1(-I)*R3(-omega)*q;
yc(4:6) = R3(-Omega)*R1(-I)*R3(-omega)*qdot;

end

function [ R ] = R1(theta)
% X-rotation
R = [1 0 0 ; 0 cos(theta) sin(theta) ; 0 -sin(theta) cos(theta)];
end


function [ R ] = R3(theta)
% Z-rotation
R = [cos(theta) sin(theta) 0 ; -sin(theta) cos(theta) 0  ; 0 0 1];
end

