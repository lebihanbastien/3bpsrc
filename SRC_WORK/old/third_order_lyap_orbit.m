%-------------------------------------------------------------------------%
% Provides a third-order approximation of a planar lyapunov orbit @ t
%-------------------------------------------------------------------------%
% @param lyap the structure containing the orbit
% @param t the initial time
% @return the initial conditions at time t in a 6*1 vector
function yv_li =  third_order_lyap_orbit(lyap, t)

%--------------------------------------------------------------------------
%Init
%--------------------------------------------------------------------------
yv_li = zeros(6,1);
gamma_i = lyap.li.gamma_i;
mu = lyap.cr3bp.mu;

%Richardson coefficients
RC = richardson_coefficients(lyap.cr3bp, lyap.li.number);

%X amplitude in Li-frame
Ax_li = lyap.Ax_estimate/gamma_i;

%Phases
Phi = 0.0; %arbitrary

%nu parameter
nu = 1+RC.s1*Ax_li^2;

%Third order estimate of the orbital period
%double Tto = (2*M_PI)/(RC.omega_p*nu);

%State estimate @ t  in the Li-frame
tau1 = RC.omega_p*nu*t + Phi;
dtau1dt = RC.omega_p*nu;

yv_li(1) = RC.a21*Ax_li^2 - Ax_li*cos(tau1) + (RC.a23*Ax_li^2)*cos(2*tau1) + (RC.a31*Ax_li^3)*cos(3*tau1);
yv_li(2) = RC.kappa*Ax_li*sin(tau1) + (RC.b21*Ax_li^2)*sin(2*tau1) + (RC.b31*Ax_li^3)*sin(3*tau1);
yv_li(3) = 0.0;

yv_li(4) = dtau1dt*( Ax_li*sin(tau1) - 2*(RC.a23*Ax_li^2)*sin(2*tau1) -3*(RC.a31*Ax_li^3)*sin(3*tau1)  );
yv_li(5) = dtau1dt*( RC.kappa*Ax_li*cos(tau1) + 2*(RC.b21*Ax_li^2)*cos(2*tau1) + 3*(RC.b31*Ax_li^3)*cos(3*tau1)   );
yv_li(6) = 0.0;

%Back into COM frame
switch(lyap.li.number) %special case of the x dimension

    case 1
        yv_li(1) = yv_li(1)*gamma_i - mu + 1 - gamma_i;

    case 2
        yv_li(1) = yv_li(1)*gamma_i - mu + 1 + gamma_i;

    case 3
        yv_li(1) = yv_li(1)*gamma_i - mu - gamma_i;

    otherwise %L2
        yv_li(1) = yv_li(1)*gamma_i - mu + 1 + gamma_i;

end

for i = 2:6 
    yv_li(i) = gamma_i*yv_li(i); %other dimensions
end
        
        
end