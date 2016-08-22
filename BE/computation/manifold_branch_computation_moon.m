function msi = manifold_branch_computation_moon(cr3bp, orbit, msi, pos, t, params, cst)
% MANIFOLD_BRANCH_COMPUTATION_MOON is equivalent to the routine 
% MANIFOLD_BRANCH_COMPUTATION_BE except that the integration is stopped at 
% the closest approach of the Moon (the smaller primary).
%
% See also MANIFOLD_BRANCH_COMPUTATION_BE.

%--------------------------------------------------------------------------
% At the closest approach to the Lunar surface, the flight path angle (fpa) is
% null. We create an 'event structure' that will help us to detect the closest
% approach as the moment when fpa = 0.0 along the trajectory.
%--------------------------------------------------------------------------
moon.event = init_event(cst.manifold.event.type.FLIGHT_PATH_ANGLE,...      %the event is triggered when the flight path angle is...
    0.0,...                                                                %equal to zero...
    cst.manifold.event.isterminal.YES,...                                  %the trajectory stops after a certain number of events
    cst.manifold.event.direction.ALL,...                                   %all direction are considered...
    cr3bp.m2.pos, cst);                                                    %the center for the computation of the fpa angle is the Moon

% Initialize the manifold with a big max_events
moon.event.max_events = 50;

% Save value of the original params structure
bool = params.plot.manifold_branch;

% Computation of the manifold branch + all the tangential points along
% the trajectory
params.plot.manifold_branch = false;
msi = manifold_branch_computation(cr3bp, orbit, msi, pos, t, params, cst, moon.event);


% Get the min distance to the surface of the Moon
mindist2moonsurface = Inf;
ind = 0;
for i = 1: size(msi.yve, 1)
    dist2moonsurface = norm(msi.yve(i,1:3) - cr3bp.m2.pos) -  cr3bp.m2.Rm/cr3bp.L;
    if(dist2moonsurface < mindist2moonsurface)
        mindist2moonsurface = dist2moonsurface;
        ind = i;
    end
end

% Recompute the manifold, stopping at the right event
params.plot.manifold_branch = bool;
moon.event.max_events = ind;
msi.ind = ind;
msi  = manifold_branch_computation(cr3bp, orbit, msi, pos, t, params, cst, moon.event);


end