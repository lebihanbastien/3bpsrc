%-------------------------------------------------------------------------%
% manifold_branch_computation(cr3bp, orbit, manifold_branch, theta, t, params, cst)
%
% Updates a manifold branch
%
% WARNING: the orbit structure should be properly computed prior to the use
% of this routine.
%-------------------------------------------------------------------------%
% Inputs:
%
% 1. cr3bp the current three-body problem
% 2. orbit the current orbit
% 3. manifold_branch the branch to be updated
% 4. theta the position of departure on the orbit. Should be taken
% in [0 1]. 
% 5. t the maximum integration time
% 6. params user parameters
% 7. cst numerical constants
%
% Outputs:
%
% 1. manifold_branch the updated manifold_branch
%
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function manifold_branch = manifold_branch_computation(cr3bp, orbit, manifold_branch, theta, t, params, cst)

%-------------------------------------------------------------------------%
%Integration until theta is reached
%-------------------------------------------------------------------------%
if(theta > 0)
    if(params.computation.type == cst.computation.MATLAB)
        %-----------------------------
        % If MATLAB routines only
        %-----------------------------
        options = odeset('Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
        [~,yvv] = ode45(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),[0 orbit.T*theta],orbit.y0,options);
        yv = yvv(end,:);
    else
        %-----------------------------
        % If MEX routines are allowed
        %-----------------------------
        [~, yv] = ode78_cr3bp(0.0, orbit.T*theta, orbit.y0, 42, cr3bp.mu);
    end
else
    yv = orbit.y0';
end

%-------------------------------------------------------------------------%
% Current STM
%-------------------------------------------------------------------------%
STM = eye(6);
for i = 1 : 6
    for j = 1 : 6
        m = 6*(i-1) + j;
        STM(i,j) = yv(m+6);
    end
end

%-------------------------------------------------------------------------%
% New vectors and initial state
%-------------------------------------------------------------------------%
if(manifold_branch.stability == cst.manifold.STABLE)
    vecp = STM*orbit.stable_direction;%New vector
    vecp_norm = vecp/norm(vecp);    
else
    vecp = STM*orbit.unstable_direction;%New vector
    vecp_norm = vecp/norm(vecp);
end

%Initial state
xs01 = (1:6)';
for i = 1 : 6
    xs01(i) = yv(end,i) + vecp_norm(1)*manifold_branch.way*cr3bp.d_man*vecp_norm(i);
end
 

%-------------------------------------------------------------------------%
%Integration
%-------------------------------------------------------------------------%
% Backwards or forward integration
if(manifold_branch.stability == 1)
    tspan = [0 -t];
else
    tspan = [0 t];
end

%-------------------------------------------------------------------------%
% Termination?
%-------------------------------------------------------------------------%
if(params.computation.type == cst.computation.MATLAB)
    %-----------------------------
    % If MATLAB routines only
    %-----------------------------
    switch(manifold_branch.event.type)
        case cst.manifold.event.type.FREE
            options = odeset('Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
            [t,yv] = ode45(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),tspan,xs01,options);
            te = t(end);
            yve = yv(end, :);
            
        case cst.manifold.event.type.ANGLE_SECTION
            options = odeset('Event',@(t,y)odezero_angle(t,y,manifold_branch.event), 'Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
            [t,yv,te,yve] = ode45(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),tspan,xs01,options);
            %If no event
            if(isempty(te))
               te = t(end);
               yve = yv(end, :);
            end
            
        otherwise
            options = odeset('Event',@(t,y)odezero_plane(t,y,manifold_branch.event), 'Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
            [t,yv,te,yve] = ode45(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),tspan,xs01,options);
            %If no event
            if(isempty(te))
               te = t(end);
               yve = yv(end, :);
            end
    end
else
    %-----------------------------
    % If MEX routines are allowed
    %-----------------------------
    switch(manifold_branch.event.type)
        case cst.manifold.event.type.FREE
            [te, yve, ~, yv] = ode78_cr3bp(0.0, tspan(2), xs01, 6, cr3bp.mu);
%         case cst.manifold.event.type.X_SECTION
%             [te, yve, ~, yv] = ode78_cr3bp(0.0, tspan(2), xs01, 6, cr3bp.mu);
        otherwise
            [te, yve, ~, yv] = ode78_cr3bp_event(0.0, tspan(2), xs01, 6, cr3bp.mu, manifold_branch.event);
            
%             %If no event
%             if(isempty(t))
%                te = t(end);
%                yve = yv(end, :);
%             end
    end
end
    

%-------------------------------------------------------------------------%
% Output
%-------------------------------------------------------------------------%
manifold_branch.termination_time = te;
manifold_branch.yv = yve;

%-------------------------------------------------------------------------%
% Compute min distance with the smallest prim
%-------------------------------------------------------------------------%
distSP = distanceToPrim(yv(1,:), orbit);
%Rm2 = cr3bp.m2.Req/cr3bp.L;
for i = 2:size(yv,1)
    if(distSP > distanceToPrim(yv(i,:), orbit))
       distSP =  distanceToPrim(yv(i,:), orbit);
    end
end

%-------------------------------------------------------------------------%
% Plotting (potentially)
%-------------------------------------------------------------------------%
if(params.plot.manifold_branch == cst.TRUE)
    manifold_plot(yv, orbit, manifold_branch, params, cst);
end


end

%-------------------------------------------------------------------------%
% Distance to the smallest primary
%-------------------------------------------------------------------------%
function distSP = distanceToPrim(yv, orbit)
    distSP =  sqrt((yv(1) - 1 + orbit.cr3bp.mu)^2 + yv(2)^2 + yv(3)^2);
end

