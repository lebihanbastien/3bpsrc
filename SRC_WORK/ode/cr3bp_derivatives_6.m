%-------------------------------------------------------------------------%
% Circular-Restricted Three-Body Problem:
%     Equations of motion for the 6-dimension state in ode45 format
%       
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function out = cr3bp_derivatives_6(t,y,mu)

%Output declaration
out = (1:6)';

%-------------------------------------------------------------------------------
% Update first derivatives of the potential \bar{U} (cf Koon et al. 2006)
%-------------------------------------------------------------------------------
d1_ub = d1_u_barre(mu,y(1),y(2),y(3));

%--------------------------------------------------------------------------
%Phase space derivatives
%--------------------------------------------------------------------------
out(1) = y(4);
out(2) = y(5);
out(3) = y(6);
out(4) = -d1_ub(1) + 2*y(5);
out(5) = -d1_ub(2) - 2*y(4);
out(6) = -d1_ub(3);

end

%-------------------------------------------------------------------------%
% First derivative of the energy potential
%-------------------------------------------------------------------------%
function xout = d1_u_barre(mu,x,y,z)
xout(1) = (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2)) - ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2)) - x;
xout(2) = (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - y - (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
xout(3) = (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
end
