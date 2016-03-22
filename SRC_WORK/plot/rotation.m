%--------------------------------------------------------------------------
% Rotation function
% Emilien Fabacher
% 2015
%--------------------------------------------------------------------------
function XYZrot = rotation(XYZ,P)
XYZrot = XYZ;
for i =1:size(XYZ,1)
    XYZrot(i,:) = (P*XYZrot(i,:)')';
end