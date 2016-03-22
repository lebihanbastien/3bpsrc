%-------------------------------------------------------------------------%
% Bicircular Four-Body Problem:
%     Equations of motion for the 6-dimension state in ode45 format
%       
% Author: BLB
% Version: 1.0
% Year: 2015
%
% @TODO: Outsource the hard coded values of ms, as, omegaS
%-------------------------------------------------------------------------%
% y = (x y z xp yp zp)
function out = bcfbp_derivatives_6(t,y,mu, omega0)

ms = 328900.54; 
as = 388.81114;
omegaS = -0.925195985520347;

%Output declaration
out = (1:6)';

%-------------------------------------------------------------------------------
% Update first derivatives of the potential \bar{U} (cf Koon et al. 2006)
%-------------------------------------------------------------------------------
d1_ub = d1_u_b4(mu, ms, as, y(1),y(2),y(3), omega0+omegaS*t);

%--------------------------------------------------------------------------
%Phase space derivatives
%--------------------------------------------------------------------------
out(1) = y(4);
out(2) = y(5);
out(3) = y(6);
out(4) = d1_ub(1) + 2*y(5);
out(5) = d1_ub(2) - 2*y(4);
out(6) = d1_ub(3);

end

%-------------------------------------------------------------------------%
% First derivative of the energy potential
%-------------------------------------------------------------------------%
function xout = d1_u_b4(mu, ms, as, x,y,z, theta)
r1 = sqrt((x+mu)^2 + y^2 + z^2);
r2 = sqrt((x-1+mu)^2 + y^2 + z^2);
rs = sqrt((x-as*cos(theta))^2 + (y-as*sin(theta))^2 + z^2);
xout(1) = x - (1-mu)/r1^3*(x+mu) - mu/r2^3*(x-1+mu) - ms/rs^3*(x-as*cos(theta)) - ms/as^2*cos(theta);
xout(2) = y - (1-mu)/r1^3*y      - mu/r2^3*y        - ms/rs^3*(y-as*sin(theta)) - ms/as^2*sin(theta);
xout(3) =   - (1-mu)/r1^3*z      - mu/r2^3*z        - ms/rs^3*z;
end