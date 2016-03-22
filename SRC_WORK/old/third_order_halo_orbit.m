function yv_li =  third_order_halo_orbit(halo, t)

%--------------------------------------------------------------------------
%Init
%--------------------------------------------------------------------------
yv_li = zeros(6,1);
gamma_i = halo.li.gamma_i;
mu = halo.cr3bp.mu;
dm = halo.dm;

%Richardson coefficients
RC = richardson_coefficients(halo.cr3bp, halo.li.number);

%X and Z-amplitude in Li-frame
Az_li = halo.Az_estimate/gamma_i;
Ax_li = sqrt(-(RC.l2*Az_li^2 + RC.Delta)/RC.l1);

%Phases
Phi = 0.0; %arbitrary

%nu parameter
nu = 1+RC.s1*Ax_li^2+RC.s2*Az_li^2;

%Third order estimate of the orbital period
%double Tto = (2*M_PI)/(RC.omega_p*nu);

%State estimate @ t  in the Li-frame
tau1 = RC.omega_p*nu*t + Phi;
dtau1dt = RC.omega_p*nu;

yv_li(1) = RC.a21*Ax_li^2 + RC.a22*Az_li^2 - Ax_li*cos(tau1) + (RC.a23*Ax_li^2 - RC.a24*Az_li^2)*cos(2*tau1) + (RC.a31*Ax_li^3 - RC.a32*Ax_li*Az_li^2)*cos(3*tau1);
yv_li(2) = RC.kappa*Ax_li*sin(tau1) + (RC.b21*Ax_li^2 - RC.b22*Az_li^2)*sin(2*tau1) + (RC.b31*Ax_li^3 - RC.b32*Ax_li*Az_li^2)*sin(3*tau1);
yv_li(3) = dm*Az_li*cos(tau1) + dm*RC.d21*Ax_li*Az_li*(cos(2*tau1)-3) + dm*(RC.d32*Az_li*Ax_li^2 - RC.d31*Az_li^3)*cos(3*tau1);

yv_li(4) = dtau1dt*( Ax_li*sin(tau1) - 2*(RC.a23*Ax_li^2 - RC.a24*Az_li^2)*sin(2*tau1) -3*(RC.a31*Ax_li^3 - RC.a32*Ax_li*Az_li^2)*sin(3*tau1)  );
yv_li(5) = dtau1dt*( RC.kappa*Ax_li*cos(tau1) + 2*(RC.b21*Ax_li^2 - RC.b22*Az_li^2)*cos(2*tau1) + 3*(RC.b31*Ax_li^3 - RC.b32*Ax_li*Az_li^2)*cos(3*tau1)   );
yv_li(6) = dtau1dt*( -dm*Az_li*sin(tau1) -2*dm*RC.d21*Ax_li*Az_li*sin(2*tau1) - 3*dm*(RC.d32*Az_li*Ax_li^2 - RC.d31*Az_li^3)*sin(3*tau1)  );


%Back into COM frame
switch(halo.li.number) %special case of the x dimension

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