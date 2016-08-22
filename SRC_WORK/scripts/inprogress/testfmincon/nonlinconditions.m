function [c,ceq] = nonlinconditions(x0, N, hL, cr3bp, msi)

c   = x0(6*N+1) - x0(7*N);
ceq = zeros(6*(N-1)+5,1);

%--------------------------------------------------------------------------
% Equality conditions
%--------------------------------------------------------------------------
for k = 1:N-1
    % State at tk
    yk   = x0((k-1)*6+1:k*6);
    tk   = x0(6*N+k);
    
    % State at tk+1
    yk1  = x0(k*6+1:(k+1)*6);
    tk1  = x0(6*N+k+1); 
    
    %Integration
    [~,pyk] = ode78_cr3bp([tk tk1], yk, cr3bp.mu);
    
    %Equality conditions
    ceq((k-1)*6+1:k*6) = pyk(end,:)' - yk1;  
end

%--------------------------------------------------------------------------
% Last 3: departure at Halo
%--------------------------------------------------------------------------
% State at t1
y1 = x0(1:6);
ceq(6*(N-1)+1:6*(N-1)+3) = y1(1:3) - msi.yv0(1:3);

%--------------------------------------------------------------------------
% Next one is: arrival at LLO
%--------------------------------------------------------------------------
% State at tN
yN = x0((N-1)*6+1:N*6);
% Equality: get on a circular orbit
ceq(6*(N-1)+4) = norm(yN(1:3)-cr3bp.m2.pos')^2 - (cr3bp.m2.Rm/cr3bp.L + hL)^2;
% Equality: tangential arrival
ceq(6*(N-1)+5) = (yN(1) - cr3bp.m2.pos(1))*(yN(4) - yN(2)) + yN(2)*(yN(5)+yN(1) - cr3bp.m2.pos(1)) + yN(3)*yN(6);

end

