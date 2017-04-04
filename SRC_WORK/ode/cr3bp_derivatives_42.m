function out = cr3bp_derivatives_42(t,y,mu)
% CR3BP_DERIVATIVES_42 provide the equations of motion for the CIRCULAR
% RESTRICTED THREE-BODY PROBLEM (CR3BP) in MATLAB ODE format.
%
% OUT = CR3BP_DERIVATIVES_42(T, Y, MU) computes the first-order
% equations of motion of the CR3BPat time T and state Y.
% The CR3BP mass ratio is MU. On top of the 6-dimensionnal state, the 36 
% equations that govern the evoluation of the State Transition Matrix (STM) 
% are also provided, for a total of 42 first-order equations of motion.
%
% The equations of motion are available in chapter 2 of Koon et al.
% "Dynamical Systems, the Three-Body Problem and Space Mission Design" 2006
% <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>. 
%
% See also CR3BP_DERIVATIVES_6
% 
% BLB 2014

%Output declaration
out = (1:42)';

%--------------------------------------------------------------------------
% Update first & second derivatives of the potential \bar{U} 
% (cf Koon et al. 2006)
%--------------------------------------------------------------------------
d1_ub = d1_u_barre(mu,y(1),y(2),y(3));
d2_ub = d2_u_barre(mu,y(1),y(2),y(3));

%--------------------------------------------------------------------------
%Phase space derivatives
%--------------------------------------------------------------------------
out(1) = y(4);
out(2) = y(5);
out(3) = y(6);
out(4) = -d1_ub(1) + 2*y(5);
out(5) = -d1_ub(2) - 2*y(4);
out(6) = -d1_ub(3);

%--------------------------------------------------------------------------
% STM derivatives
%--------------------------------------------------------------------------
%STM is updated
Phi = vectorToMatrix(y, 6, 6, 6);

% Omega 3*3
Omega = zeros(3);
Omega(1,2) =  1;
Omega(2,1) = -1;

%Super matrix Df(x) (6*6)
Df = [zeros(3) eye(3) ; -d2_ub 2*Omega];

%Derivative
d_Phi = Df * Phi;

%Storage in output
out = matrixToVector(out, d_Phi, 6, 6, 6);
% for i = 1 : 6
%     for j = 1 : 6
%         m = 6*(i-1) + j;
%         out(6+m) = d_Phi(i,j);
%     end
% end

end

%-------------------------------------------------------------------------%
% First derivative of the energy potential
%-------------------------------------------------------------------------%
function xout = d1_u_barre(mu,x,y,z)

r1 = ((mu + x - 1)^2 + y^2 + z^2);
r2 = ((mu + x)^2 + y^2 + z^2);

xout(1) = mu*(x - 1 + mu)/r1^(3/2) + (1 - mu)*(x + mu)/r2^(3/2) - x;
xout(2) = mu*y/r1^(3/2)            + (1 - mu)*y/r2^(3/2)        - y;
xout(3) = mu*z/r1^(3/2)            + (1 - mu)*z/r2^(3/2);

end

%-------------------------------------------------------------------------%
% Second derivative of the energy potential
%-------------------------------------------------------------------------%
function xout = d2_u_barre(mu,x,y,z)

r1 = ((x - 1 + mu)^2 + y^2 + z^2);
r2 = ((x + mu)^2 + y^2 + z^2);


xout(1,1) =  mu/r1^(3/2) + (1 - mu)/r2^(3/2) - 3*mu*(x - 1 + mu)^2/r1^(5/2) - 3*(x + mu)^2*(1 - mu)/r2^(5/2) - 1;
xout(1,2) =  3*y*(mu + x)*(mu - 1)/r2^(5/2)  - 3*mu*y*(x - 1 + mu)/r1^(5/2);
xout(1,3) =  3*z*(mu + x)*(mu - 1)/r2^(5/2)  - 3*mu*z*(x - 1 + mu)/r1^(5/2);

xout(2,1) = xout(1,2);
xout(2,2) = mu/r1^(3/2) - (mu - 1)/r2^(3/2) + 3*y^2*(mu - 1)/r2^(5/2) - 3*mu*y^2/r1^(5/2) - 1;
xout(2,3) = 3*y*z*(mu - 1)/r2^(5/2) - 3*mu*y*z/r1^(5/2);

xout(3,1) = xout(1,3);
xout(3,2) = xout(2,3);
xout(3,3) = mu/r1^(3/2) - (mu - 1)/r2^(3/2) + 3*z^2*(mu - 1)/r2^(5/2) - 3*mu*z^2/r1^(5/2);

end

