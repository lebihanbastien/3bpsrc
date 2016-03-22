%--------------------------------------------------------------------------
% Example n°5: this matlab file makes use of the richardson third order
% approximation to build an EML2 planar lyapunov orbit.
% Then, its exterior stable and unstable manifolds are computed
% and integrated up to given termination conditions:
% - The stable manifold is stopped at a given angle with respect to the
% Earth-Moon line, centered at the Earth.
% - The unstable manifold is stopped at x = x0.
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% User input data
user.Az      = 8000;             %size of the orbit in [km]
user.sangle  = degtorad(45);     %termination angle for the stable manifold (see events) in [rad]
user.x0      = 0.0;              %termination value for the unstable manifold (see events) in [adim]

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Orbit
orbit = init_lyap_orbit(cr3bp, cr3bp.l2, user.Az, cst);

%% Orbit computation
orbit = lyapunov_orbit_computation(cr3bp, orbit, default, cst);



%% Events
stable.event = init_event(cst.manifold.event.type.ANGLE_SECTION,...        %the event is triggered at a given angle...
                          user.sangle,...                                  %given in user data...
                          cst.manifold.event.isterminal.YES,...            %the trajectory stops at the first ocurrence...
                          cst.manifold.event.direction.ALL,...             %all direction are considered...
                          cr3bp.m1.pos, cst);                              %the center for the computation of the angle is the Earth
               
unstable.event = init_event(cst.manifold.event.type.X_SECTION,...          %the event is triggered at a given value x = x0
                            user.x0,...                                    %given in user data...
                            cst.manifold.event.isterminal.YES,...          %the trajectory stops at the first ocurrence...
                            cst.manifold.event.direction.ALL,...           %all direction are considered...
                            cr3bp.m1.pos, cst);                            %the center for the computation of the angle is the Earth               

%% Manifold initialization
% Init
                                      
% The stable manifold is stopped @45° wrt the Earth
manifold_branch_stable = init_manifold_branch_event(cst.manifold.STABLE, ...
                                                    cst.manifold.EXTERIOR,...
                                                    stable.event);
                                                
% The stable manifold is stopped @x = 0
manifold_branch_unstable = init_manifold_branch_event(cst.manifold.UNSTABLE, ...
                                                      cst.manifold.EXTERIOR,...
                                                      unstable.event);

%% Manifold computation
t = 20;

% Stable
for theta = 0:0.1:1 %position on the orbit
    manifold_branch_stable = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable, theta, t, default, cst);
end

% Unstable
for theta = 0:0.1:1 %position on the orbit
    manifold_branch_unstable = manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable, theta, t, default, cst);
end