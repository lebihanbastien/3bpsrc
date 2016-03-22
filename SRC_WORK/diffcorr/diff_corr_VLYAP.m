%--------------------------------------------------------------------------
%Differential correction for 3D periodic orbit computation
%Suitable for e.g. halo orbits
% - cst.corr.Z0_FIXED has to be used for small orbits
% - cst.corr.X0_FIXED has to be used for bigger orbits
%--------------------------------------------------------------------------
% @param yv0 the approximated initial conditions
% @param cr3bp the structure containing the parent CR3BP
% @param orbit the structure containing the parent orbit
% @param params the structure containing the comptuation parameters
% @param cst the structure containing the numerical constants
% @return the udpated structure orbit
function orbit = diff_corr_VLYAP(yv0, cr3bp , orbit, params, cst)
%Iteration counts
iter = 0;

%Event structure
val_par = init_event(cst.manifold.event.type.Y_SECTION,...
                     0.0,...
                     cst.manifold.event.isterminal.YES,...
                     cst.manifold.event.direction.ALL,... 
                     cr3bp.m1.pos,...
                     cst);
                 
%Options for Matlab integration
options = odeset('Events',@odezero_y,'Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);          

%Diff corr loop
while(true)
    iter = iter+1;
    
    %----------------------------------------------------------------------
    % Integration stops @y=0
    %----------------------------------------------------------------------
    if(params.computation.type == cst.computation.MATLAB)
        %-----------------------------
        % If MATLAB routines only
        %-----------------------------
        [~,yv,te,ve,~] = ode45(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),[0 10],yv0,options);
    else
        %-----------------------------
        % If MEX routines are allowed
        %-----------------------------
        if(params.plot.diff_corr == 1)
            [te, ve, ~, yv] = ode78_cr3bp_event(0.0, 10, yv0, 42, cr3bp.mu, val_par);
        else
            [te, ve] = ode78_cr3bp_event(0.0, 10, yv0, 42, cr3bp.mu, val_par);
        end
        
    end

    
    %----------------------------------------------------------------------
    % Plotting (potentially)
    %----------------------------------------------------------------------
    if(params.plot.diff_corr == 1) %plotting
        figure(3)
        if(iter==1)
            plot(yv(:,2)*cr3bp.L*10^(-3),yv(:,3)*cr3bp.L*10^(-3), 'g');
        else
            plot(yv(:,2)*cr3bp.L*10^(-3),yv(:,3)*cr3bp.L*10^(-3), 'r');
        end
        hold on
        xlabel('Y (x 10^3 km)')
        ylabel('Z (x 10^3 km)')
        axis equal
        grid on
        title('Halo orbit in Y-Z plane');
        
    elseif(params.plot.diff_corr==2 && iter==1) %only first iter
        
        figure(3)
        plot(yv(:,2)*cr3bp.L*10^(-3),yv(:,3)*cr3bp.L*10^(-3), 'r');
        hold on
        xlabel('Y (x 10^3 km)')
        ylabel('Z (x 10^3 km)')
        axis equal
        grid on
        title('Halo orbit in Y-Z plane');      
    end
    
    %----------------------------------------------------------------------
    % Correction
    %----------------------------------------------------------------------
    
    % If corr_type == cst.corr.Z0_FIXED:
    %                      Af                         ppf           Bf
    % [dzf]       = ( [Phif31  Phif35] - 1/ypoint * [zp]  * [Phi21   Phi25] ) * [   dx0  ]
    % [dxfpoint]      [Phif41  Phif45]              [xpp]                       [dy0point]
    %
    
    % If corr_type == cst.corr.X0_FIXED:
    %                      Af                         ppf           Bf
    % [dzf]       = ( [Phif33  Phif35] - 1/ypoint * [zp]  * [Phi23   Phi25] ) * [   dz0  ]
    % [dxfpoint]      [Phif43  Phif45]              [xpp]                       [dy0point]
    %
    
    %Af and Bf
    Af = zeros(2);
    Bf = (1:2);
    if(params.diff_corr.type == cst.corr.Z0_FIXED)
        Af(1,1) = ve(6+13); %Phif31
        Af(1,2) = ve(6+17); %Phif35
        Af(2,1) = ve(6+19); %Phif41
        Af(2,2) = ve(6+23); %Phif45
        
        Bf(1) = ve(6+7);    %Phif21
        Bf(2) = ve(6+11);   %Phif25
    else
        Af(1,1) = ve(6+15); %Phif33
        Af(1,2) = ve(6+17); %Phif35
        Af(2,1) = ve(6+21); %Phif43
        Af(2,2) = ve(6+23); %Phif45
        
        Bf(1) = ve(6+9);    %Phif23
        Bf(2) = ve(6+11);   %Phif25
    end
    
    
    %State subvector
    y_state = (1:6)';
    for i = 1 : 6
        y_state(i) = ve(i);
    end
    
    %Derivative of y_state @ y=0 in y_state_p
    y_state_p = cr3bp_derivatives_6(te, y_state, cr3bp.mu);
       
    %ppf vector
    ppf = (1:2)';
    ppf(1) = y_state_p(3)/y_state(5);  %zp
    ppf(2) = y_state_p(4)/y_state(5);  %xfpp
    
    
    %New Af
    Af = Af - ppf*Bf;
     
    %Inversion of the error
    dyf = [-y_state(3) ; -y_state(4)];  % = [- zf ; -xfpoint]
    
    %Stops if precision is good enough
    if(abs(dyf(1)) < params.diff_corr.precision &&  abs(dyf(2)) < params.diff_corr.precision);
        orbit.T12 = te; %1/2 period
        orbit.T = 2*te; %period
        break;
    end
    
    if(iter > 50)
        disp('WARNING: maximum iterations reached in differential_correction');
        break;
    end
    
    %Sinon, calcul de dx0
    dy0 = Af \ dyf;
    
    %Updating initial state
    if(params.diff_corr.type == cst.corr.Z0_FIXED)
        yv0(1) = yv0(1) + dy0(1);
        yv0(5) = yv0(5) + dy0(2);  
    else
        yv0(3) = yv0(3) + dy0(1);
        yv0(5) = yv0(5) + dy0(2);
    end
   
end

%Orbit update
orbit.y0 = yv0;
end

