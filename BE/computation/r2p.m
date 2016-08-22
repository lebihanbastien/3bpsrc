function [ p ] = r2p( r, I, Omega )
    p = R1(I)*R3(Omega)*r;
end

function [ R ] = R1(theta)
% X-rotation
R = [1 0 0 ; 0 cos(theta) sin(theta) ; 0 -sin(theta) cos(theta)];
end


function [ R ] = R3(theta)
% Z-rotation
R = [cos(theta) sin(theta) 0 ; -sin(theta) cos(theta) 0  ; 0 0 1];
end

