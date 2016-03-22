%-------------------------------------------------------------------------%
% Circular-Restricted Three-Body Problem:
%     Equations of motion for the 6-dimension state in ode45 format
%       
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function out = cr3bp_mn_6(t,y,mu,c1,gamma)

%Output declaration
out = (1:6)';

%--------------------------------------------------------------------------
%Phase space derivatives
%--------------------------------------------------------------------------
out(1) = y(4) + y(2);
out(2) = y(5) - y(1);
out(3) = y(6);

e = c1 - mu/gamma;
m = c1 - (mu-1)/gamma;
qpe2 = (y(1) - e)^2 + y(2)^2 + y(3)^2;
qpm2 = (y(1) - m)^2 + y(2)^2 + y(3)^2;


out(4) = +y(5)-c1...
         -1/gamma^3*( (1-mu)/qpe2^(3/2)*(y(1) - e) + mu/qpm2^(3/2)*(y(1) - m) );
     
out(5) = -y(4)...
         -1/gamma^3*( (1-mu)/qpe2^(3/2)*y(2) + mu/qpm2^(3/2)*y(2) );
     
out(6) = -1/gamma^3*( (1-mu)/qpe2^(3/2)*y(3) + mu/qpm2^(3/2)*y(3) );
     
    

end
