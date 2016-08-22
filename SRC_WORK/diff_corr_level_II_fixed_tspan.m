function [tmdn, ymdn, yma] = diff_corr_level_II_fixed_tspan(tmd, ymd, N2, dcs, inputType)

%--------------------------------------------------------------------------
% Copy the departure state
%--------------------------------------------------------------------------
tmdn = tmd;
ymdn = ymd;

%--------------------------------------------------------------------------
% Loop
%--------------------------------------------------------------------------
normC = 1e5;
iter = 0;
while(normC > 1e-6 && iter < 5)
    
    fprintf('---------------------------------------------------\n');
    fprintf('Level I\n');
    fprintf('---------------------------------------------------\n');
    
    %----------------------------------------------------------------------
    % Level one. After this point: continuity in position is ensured
    %----------------------------------------------------------------------
    [tmdn, ymdn, yma] = diff_corr_level_I_fixed_tspan(tmdn, ymdn, N2, dcs, inputType);
    
    fprintf('---------------------------------------------------\n');
    fprintf('Level II\n');
    fprintf('---------------------------------------------------\n');
    
    
    %----------------------------------------------------------------------
    % Initialize the super matrix & the error vector
    %----------------------------------------------------------------------
    % Number of inner points
    ni = N2-2;
    % Initialization
    Ms = zeros(3*ni, 12*ni + 4*(ni-1));
    dv = zeros(3*ni, 1);
    
    %----------------------------------------------------------------------
    % Level two
    %----------------------------------------------------------------------
    for k = 1:N2-2 %only the inner points!
        
        %------------------------------------------------------------------
        % The arrival state contains all the desired elements
        %------------------------------------------------------------------
        % Final state
        ypa = yma(k+1,:)';
        tpa = tmdn(k+1);
        % STM
        STMa = vectorToMatrix(ypa, 6, 6, 6);
        % Derivatives at tpa
        ypadot = qbcp_vfn_novar(tpa, ypa(1:6), dcs)';
        
        %------------------------------------------------------------------
        % The departure state contains all the desired elements
        %------------------------------------------------------------------
        % Final state
        ypd = ymdn(k+1,:)';
        tpd = tmdn(k+1);
        % STM
        STMd = inv(vectorToMatrix(yma(k+2, :)', 6, 6, 6));
        % Derivatives at tf
        ypddot = qbcp_vfn_novar(tpd, ypd(1:6), dcs)';
        
        %------------------------------------------------------------------
        % Level II matrices
        %------------------------------------------------------------------
        M0  = STMa(4:6,4:6)*inv(STMa(1:3,4:6))*STMa(1:3, 1:3) - STMa(4:6,1:3);
        Mt0 = ypadot(4:6) - STMa(4:6,4:6)*inv(STMa(1:3,4:6))*ypadot(1:3);
        Mp  = STMd(4:6,4:6)*inv(STMd(1:3,4:6)) - STMa(4:6,4:6)*inv(STMa(1:3,4:6));
        Mtp = STMa(4:6,4:6)*inv(STMa(1:3,4:6))*ypadot(1:3) - ypadot(4:6) - STMd(4:6,4:6)*inv(STMd(1:3,4:6))*ypddot(1:3) + ypddot(4:6);
        Mf  = STMd(4:6,1:3) - STMd(4:6,4:6)*inv(STMd(1:3,4:6))*STMd(1:3, 1:3);
        Mtf = STMd(4:6,4:6)*inv(STMd(1:3,4:6))*ypddot(1:3) - ypddot(4:6);
        
        
        %------------------------------------------------------------------
        % Update the super matrix
        %------------------------------------------------------------------
        kl = 3*(k-1);
        kc = 4*(k-1);
        
        Ms(1+kl:3+kl, 1+kc:3+kc) = M0;
        Ms(1+kl:3+kl, 4+kc) = Mt0;
        
        Ms(1+kl:3+kl, 5+kc:7+kc) = Mp;
        Ms(1+kl:3+kl, 8+kc) = Mtp;
        
        Ms(1+kl:3+kl, 9+kc:11+kc) = Mf;
        Ms(1+kl:3+kl, 12+kc) = Mtf;
        
        %------------------------------------------------------------------
        % Update the error vector
        %------------------------------------------------------------------
        dv(1+kl:3+kl) = -(ypd(4:6) - ypa(4:6));
    end
    
    %----------------------------------------------------------------------
    % Inverse the level II system
    %----------------------------------------------------------------------
    dr = Ms'*(Ms*Ms'\dv);
    
    %----------------------------------------------------------------------
    % Norm
    %----------------------------------------------------------------------
    normC = norm(dv);
    fprintf('normC = %5.15e\n', normC);
    fprintf('norm(dr) = %5.15e\n', norm(dr));
    
    %----------------------------------------------------------------------
    % Update the state
    %----------------------------------------------------------------------
    for k = 1:N2
        kc = 4*(k-1);
        if(k ~= 1 && k ~= N2) 
            ymdn(k,1:3) = ymdn(k, 1:3) + dr(1+kc:3+kc)';
            tmdn(k)     = tmdn(k) + dr(4+kc);
        end
    end
    
    %----------------------------------------------------------------------
    % Update number of iterations
    %----------------------------------------------------------------------
    iter = iter + 1;
    
end

end
