%--------------------------------------------------------------------------
% jacobi(yv, mu). Jacobi constant computation.
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------
function C = jacobi(yv, mu)
if(size(yv,2) == 3) %if only position is given, velocity is supposed null
    C = -2*u_barre(yv, mu);
else
    C = -2*u_barre(yv, mu) - (yv(4)^2 + yv(5)^2 +yv(6)^2);
end
end

%--------------------------------------------------------------------------
% Energy potential \bar{U} (see Koon et al. 2008)
%--------------------------------------------------------------------------
function xout = u_barre(yv, mu)
mu1 = 1 - mu;
mu2 = mu;
r1 = sqrt( (yv(1)+mu2)^2 + yv(2)^2 + yv(3)^2 );
r2 = sqrt( (yv(1)-mu1)^2 + yv(2)^2 + yv(3)^2 );
xout = -1/2*(yv(1)^2+yv(2)^2) - mu1/r1 - mu2/r2 - 1/2*mu1*mu2;
end