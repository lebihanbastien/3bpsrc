%--------------------------------------------------------------------------
% Advanced example n°2:
%
% This matlab file compute an Earth-to-Halo transfer with:
% - A first maneuver @LEO
% - A second maneuver at a given position on the manifold, when the flight
% path angle wrt the Earth is null (perigee/apogee-point scenario)
%
% The transfer is computed for several points
% on the halo orbit and, in the corresponding loop, the previous DeltaV
% is used to compute the transfer for the next point,
% providing that the previous solution is valid
%
% Both maneuvers are computed backwards in time, starting from a given
% point on an EML1,2-halo orbit. The LEO is defined only by its altitude
% (the final output.earth.inclination is governed by the correction scheme)
%
% The algorithm of the differential correction prodedure @LEO has been
% inpired by Gordon 2008 (Master thesis)
%
% This matlab file makes of the abacus
% ./data/halo_init_matrix_EML1,2.dat to generate an EML1,2 halo orbit and its
% exterior stable and unstable manifolds.
%
% WARNING: this computation includes:
% 1. A procedure to leave the close vicinity of the Moon
% 2. A first guess fo the deltaV of the second maneuver
%    (at the insertion point of the stable/unstable manifold branch)
% Both are arbitrarily implemented from heuristic considerations. They
% both need to be improved since the final result greatly depends on them.
%
% Possible changes:
% -1. Only the solutions with eastward LEO rotation should be
%    kept.
% 0. Select the events that are far enough from the initial orbit to avoid
% false perigee/apogee detection.
% 1. Look at Alessi et al. 2009 to get a better differential corrector?
% For now, it is the same as the lunar flyby case, which may be not
% suitable in this more general case.
% 2. Use also hohmann like transfer to compute a first guess. Should
% improve apogee results which are very bad right now.
% 3. Use true lambert problem resolution and make the tof + LEO position (angle) vary.
% 4. Using decreasing/increasing and the maximum number of events should
% allow to get a precise apogee/perigee number. For now, the user has to do
% it manually.
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters
% Integration precision values are dramatically lowered to allow fast
% computation: only for plotting!
default.plot.XY = true;
default.plot.firstPrimDisp = cst.TRUE;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
%Orbit
orbit = init_orbit(cr3bp, cr3bp.l1,  cst.orbit.type.HALO, cst.orbit.family.NORTHERN, 15clos000, cst);
%Interpolation matrix
halo_init = halo_init_EML1;

%% User input data
user.hLEO       = 185;           %Desired LEO altitude [km]
user.theta      = 0;             %Arbitrary position on the orbit, in [0 1]
user.showSteps  = false;         %To show/hide the steps of the corrective scheme

user.flyby.tangentManeuver = false;     %Is the Flyby Maneuver forced to be tangent? Very restrictive! May lead to no solution
user.flyby.fbangle  = degtorad(-130);   %Flyby angle at the moon [rad]


%Computed from user inputs, or arbitrary
user.hLEOa = user.hLEO/cr3bp.L;  %Desired LEO altitude [adim]
user.t0    = 3.5;                 %Integration duration (arbitrary, will not be reached because of the termination @lunar flyby


%% Orbit computation
orbit = halo_orbit_interpolation(cr3bp, orbit, halo_init, default, cst);
%orbit = orbit_computation(cr3bp, orbit, default, cst);

%% Events
earth.event = init_event(cst.manifold.event.type.FLIGHT_PATH_ANGLE,...     %the event is triggered when the flight path angle is...
                         0.0,...                                           %equal to zero...
                         cst.manifold.event.isterminal.NO,...             %the trajectory stops at the first ocurrence...
                         cst.manifold.event.direction.INCREASING,...      %increasing direction are considered (perigee)...
                         cr3bp.m1.pos, cst);                               %the center for the computation of the angle is the Moon

% More than one event are considered. We use the last one to compute the
% bridge leg to LEO
earth.event.max_events = 2;

%% Additional options
%Associated ode options
earth.options = odeset('Events',@(t,y)odezero_flightpathangle(t,y,earth.event),...
    'Reltol', default.ode45.RelTol,...
    'Abstol', default.ode45.AbsTol);

% Position of the primaries
moon.position  = cr3bp.m2.pos;
earth.position = cr3bp.m1.pos;

%% Manifold initialization

% Event: fpa = 0
% manifold_branch_stable  = init_manifold_branch_event(cst.manifold.STABLE,...    %the manifold is taken stable...
%                                                      cst.manifold.INTERIOR,...  %and interior...
%                                                      earth.event);              %the integration is stopped when a certain angle wrt to the moon is reached
                                                 
% No event                                                
manifold_branch_stable  = init_manifold_branch_event(cst.manifold.STABLE,...
                                                     cst.manifold.INTERIOR);
%% Loop on the position in L2
%Waitbar
h = waitbar(0,'Computation in progress...');
%--------------------------------------------------------------------------
% A rough first guess is used in lfb for the first entry in the loop.
% Then, output.flyby.deltaV can be used as first guess for the next
% iteration.
%--------------------------------------------------------------------------
isPreviousSolution = false;

% For saving outputs throughout the loop
it = 1;

%Change earth event to match needs on earth2manifold arc
earth.event.direction  = cst.manifold.event.direction.ALL;
earth.event.max_events = 1;

for theta = 0:0.01:1
    %----------------------------------------------------------------------
    %New starting point
    %----------------------------------------------------------------------
    user.theta = theta;
    %----------------------------------------------------------------------
    %Manifold computation
    %----------------------------------------------------------------------
    manifold_branch_stable = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable, user.theta, user.t0, default, cst);
    
    %----------------------------------------------------------------------
    % Lambert arc
    %----------------------------------------------------------------------
    if(isPreviousSolution) % a first guess is added from a previous solution
        [output, isPreviousSolution] = earth2manifold(manifold_branch_stable, cr3bp, earth, moon, user, default, cst, output.flyby.deltaV);
    else % no previous solution
        [output, isPreviousSolution] = earth2manifold(manifold_branch_stable, cr3bp, earth, moon, user, default, cst);
    end
    
    %----------------------------------------------------------------------
    % Saved outputs throughout the loop
    %----------------------------------------------------------------------
    if(isPreviousSolution)
        deltaVD(it) = output.deltaV_dim;
        thetaV(it)   = user.theta;
        Ttot(it)    = output.Ttot_dim;
        it = it+1;
    end
    
    %----------------------------------------------------------------------
    % Waitbar
    %----------------------------------------------------------------------
    waitbar(theta);
end

close(h)

%% Maneuver cost vs position on the Halo orbit
figure;
hold on
grid on
xlabel('Position on the Halo orbit in [0 1]');
ylabel('Total maneuver cost [km/s]');
plot(thetaV, deltaVD, 'ob', 'MarkerFaceColor', 'black');


%% Maneuver cost vs time of flight
figure;
hold on
grid on
xlabel('Time of flight [days]');
ylabel('Total maneuver cost [km/s]');
plot(Ttot, deltaVD, 'ob', 'MarkerFaceColor', 'black');

