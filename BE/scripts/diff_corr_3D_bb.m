function orbit = diff_corr_3D_bb(y0, cr3bp , orbit, params, cst)
% DIFF_CORR_3D_BB differential correction to compute 3D periodic orbits of
% the CRTBP, symmetric with respect to the xz-plane (halo orbits only, in
% this simplified case).
%
% DIFF_CORR_3D_BB(Y0, CR3BP, ORBIT, PARAMS, CST) iteratively correct the 
% initial state Y0 via an iterative Newton method. This method is 
% illustrated below.
%
% For Halo orbits, at each step:
% The free variables are
%       X0 = [x0 vy0]^T
%
% The constraint vector is FX = F(X) = [y vx vz]^T. The goal is to
% satisfy the constraint FX = 0, so as to produce a periodic orbit
% symmetric with respect to the xz-plane, by correcting X0.
%
% Starting with initial condition Y0 = [x0 y0 z0 vx0 vy0 vz0]^T,
% the equations of motion of the CRTBP are integrated until the xz-plane is
% reached. With this choice the first constraint y = 0 is automatically 
% satisfied.
%
% Then, the Newton's correction dX0 to be applied so as to obtain
%               FXr = [vx vz]^T = 0
% satisfies the following equation:
%
%   FXr                  Af               ppf               Bf        dX0
% [dvx]  = ( [Phi41  Phi45] - 1/ypoint * [ax] * [Phi21   Phi25] ) * [dx0 ]
% [dvz]      [Phi61  Phi65]              [az]                       [dvy0]
%
%
% where Phi is the State Transition Matrix (STM), numerically
% integrated along with the state.
%
% This correction is applied iteratively on the initial state until a
% use-defined precision threshold is reached.
%
% See Howell 1984, for details <a href="matlab:
% web('http://adsabs.harvard.edu/full/1984CeMec..32...53H','-browser')">(link)</a>.
%
%
% At the end of the routine, the following elements are updated:
% - orbit.y0:         initial conditions.
% - orbit.T12:        half period.
% - orbit.T:          period.
%
% BLB 2016.

%--------------------------------------------------------------------------
%Iteration counts
%--------------------------------------------------------------------------
iter = 0;

%--------------------------------------------------------------------------
%Type of event: reaching the y = 0 plane
%--------------------------------------------------------------------------
%Associated event structure for MEX routines
val_par = init_event(cst.manifold.event.type.Y_SECTION,...
                     0.0,...
                     cst.manifold.event.isterminal.YES,...
                     cst.manifold.event.direction.ALL,...
                     cr3bp.m1.pos,...
                     cst);

%--------------------------------------------------------------------------
% Differential correction loop
%--------------------------------------------------------------------------
while(true)
    iter = iter+1;
    
    if(iter > 50)
        disp('WARNING: maximum iterations reached in differential_correction');
        break;
    end
    
    %----------------------------------------------------------------------
    % Integration from y0, stops at y = 0:
    %   - te  is the final time
    %   - yve is the final 42-dimensionnal state
    %   - yv  is the entire trajectory
    %----------------------------------------------------------------------
    [te, yve, ~, yv] = ode78_cr3bp_event([0 10], y0, cr3bp.mu, val_par);

    %----------------------------------------------------------------------
    % Update the final state: FX = [vx vz]^T.
    %----------------------------------------------------------------------
    %
    % A COMPLETER (~ 1 Ligne)
    %
    
    %----------------------------------------------------------------------
    % Compute the first order correction with Newton's method
    %----------------------------------------------------------------------
    %  FXr               Af                   ppf         Bf              dX0
    % [dvx]  = ( [Phi41  Phi45] - 1/ypoint * [ax] * [Phi21   Phi25] ) * [dx0 ]
    % [dvz]      [Phi61  Phi65]              [az]                       [dvy0]
    %----------------------------------------------------------------------
    %
    % A COMPLETER (~ 20 Lignes)
    %
    % Aide : la fonction cr3bp_derivatives_6 fournit la dérivée temporelle
    % de l'état yv = [x y z vx vy vz]^T.
    %
    
    %----------------------------------------------------------------------
    % The first order correction is computed
    %----------------------------------------------------------------------
    %
    % A COMPLETER (~ 1 Ligne)
    %
    
    %----------------------------------------------------------------------
    %Updating initial state
    %----------------------------------------------------------------------
    %
    % A COMPLETER (~ 2 Lignes)
    %
    
    
    %----------------------------------------------------------------------
    % Plotting (potentially)
    %----------------------------------------------------------------------
    if(params.plot.diff_corr)
        plotDiffCorr(yv, cr3bp, iter);
    end
    
end

%Orbit update
orbit.y0  = y0;   %initial state (corrected)
orbit.T12 = te;   %1/2 period
orbit.T   = 2*te; %period

end


%--------------------------------------------------------------------------
% Plotting the iterations
%--------------------------------------------------------------------------
function [] = plotDiffCorr(yv, cr3bp, iter)
figure(1)
hold on
grid on
xlabel(strcat('X', ' (km)'));
ylabel(strcat('Y', ' (km)'));
title('Differential correction process');
if(iter==1)
    plot(yv(:,1)*cr3bp.L,yv(:,2)*cr3bp.L, 'g');
else
    plot(yv(:,1)*cr3bp.L,yv(:,2)*cr3bp.L, 'r');
end
end