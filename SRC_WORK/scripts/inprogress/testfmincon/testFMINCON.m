%--------------------------------------------------------------------------
% This matlab file allows to build:
% - an EML2 planar lyapunov orbit (~8000 km of radius)
% - an EML2 halo orbit (~10 000km of vertical extension)
% - an EML2 vertical lyapunov orbit (~30 000km of vertical extension)
%
% BLB 2016
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Same for a halo orbit
%Initialization
halo = init_orbit(cr3bp, ...       % Parent CR3BP
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.HALO, ...       % HALO orbit
    cst.orbit.family.NORTHERN, ... % Northern class
    10000, ...                     % Of vertical extension ~ 10000 km
    cst);                          % Numerical constants

%Computation
halo = orbit_computation(cr3bp, halo, default, cst);

%% Integration duration: arbitrarily fixed to 22 days
t = 22*cst.env.days*2*pi/cr3bp.T;

%% Initialization of the manifold
msi = init_manifold_branch(cst.manifold.UNSTABLE,...    %it is an UNSTABLE manifold.
                           cst.manifold.INTERIOR);      %it is an INTERIOR manifold: it will go towards the Moon.
                       
%% Computation of the manifold
% Departure position on the orbit in the interval [0,1]
theta = 0.0;
% Computation of the manifold branch
msi = manifold_branch_computation_moon(cr3bp,...   %parent CRTBP
                                       halo,...   %parent orbit
                                       msi,...     %current manifold branch
                                       theta,...   %departure position on the orbit
                                       t,...       %time of flight on the manifold branch
                                       default,... %default parameters
                                       cst);       %constants
                                   
%% Preparation to FMINCON

% 1. Get N points along the trajectory
N = 60;
tspan = linspace(0.0, msi.te(end), N);

% 2. Integrate to get the state on these N points
options = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol);
[tvv,yvv] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),tspan,msi.yv0,options);

%3. Plot them on Figure 4
figure(4)
hold on
%plot3(yvv(:,1)*cr3bp.L, yvv(:,2)*cr3bp.L, yvv(:,3)*cr3bp.L, '*');



%% Fmincon
%options = optimoptions('fmincon','Display','iter','Algorithm', 'interior-point', 'MaxFunEvals', 5000);
options = optimoptions('fmincon','Display','iter','Algorithm', 'sqp', 'MaxFunEvals', 10000);

% Build the state x0 = [y1,...,yN, t1,...tN]
x0 = [];
for i = 1:N
    x0 = [x0 ; yvv(i,:)'];
end
for i = 1:N
    x0 = [x0 ; tvv(i)];
end

%% Initial altitude is?
hL = norm(yvv(end, 1:3) - cr3bp.m2.pos) - cr3bp.m2.Rm/cr3bp.L - 1e-2;

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

nonlcon = @(x)nonlinconditions(x, N, hL, cr3bp, msi);
fun     = @(x)minconditions(x, N, msi, cr3bp, hL);
x       = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

for k = 1:N
    % State at tk
    yvv2(k,:) = x((k-1)*6+1:k*6)';
end

figure(4)
hold on
plot3(yvv2(:,1)*cr3bp.L, yvv2(:,2)*cr3bp.L, yvv2(:,3)*cr3bp.L, 'o');

%% Replot
for k = 1:N-1
    [~,~, ~, yv2] = ode78_cr3bp([x(6*N+k) x(6*N+k+1)], yvv2(k,:), cr3bp.mu);

    figure(4)
    hold on
    plot3(yv2(:,1)*cr3bp.L, yv2(:,2)*cr3bp.L, yv2(:,3)*cr3bp.L);

end

%% Maneuver at Halo
DV = norm(msi.yv0(4:6) - yvv2(1,4:6)')*cr3bp.L*2*pi/(cr3bp.T)*1000;
figure(4)
hold on
title(['DVh = ', num2str(DV, '%5.5f'), ' m/s']);

%% Plot orbit to target
vl = cr3bp.m2.pos*cr3bp.L;
rl = cr3bp.m2.Rm + hL*cr3bp.L;

figure(4);
hold on
[xS, y, z] = sphere(16);
h = surf(rl * xS + vl(1), rl * y, rl * z);
set(h, 'FaceAlpha', 0.1, 'EdgeColor', 0.5*ones(3,1), 'FaceColor', 0.5*ones(3,1));
shading interp