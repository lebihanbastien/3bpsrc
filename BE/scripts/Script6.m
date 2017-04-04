% Sixth script of the BE.
% Computing a Halo orbit around EML2 + its unstable directions + a branch
% of the unstable manifold, stopping at the closest approach of the lunar
% surface + targeting a specific altitude with a Hohmann transfer
%
% 2016

%% Call of the first script: building the orbit
Script1

%% Stop plotting on figure 1
default.plot.XY = false;

%% Integration duration of the manifold arbitrarily fixed to 5 days
t = 22*cst.env.days*2*pi/cr3bp.T;

%% Initialization & computation of the manifold
msi = init_manifold_branch(cst.manifold.UNSTABLE, cst.manifold.INTERIOR);
% Departure position on the orbit in [0, 1]
theta = 0.2; %0.23 for 50 km
% Computation of the manifold branch.
msi = manifold_branch_computation_moon(cr3bp, orbit, msi, theta, t, default, cst);

%% Corresponding LLO

%--------------------------------------------------------------------------
% Initialize some useful constants
%--------------------------------------------------------------------------
% Lunar radius
rm = cr3bp.m2.Rm/cr3bp.L;
% Lunar position
pm = cr3bp.m2.pos;
% Earth-Moon distance
Lf = cr3bp.L;
% The velocity factor
Tf = cr3bp.L*2*pi/cr3bp.T;
% Maximum altitude allowed at the Moon : 500 km.
maxalt = 500/cr3bp.L;

%--------------------------------------------------------------------------
% Current state
%--------------------------------------------------------------------------
yv = msi.yve(end, :)';
tv = msi.te(end);

% Current state wrt the Moon
yvm = yv - [pm 0 0 0]';

% Altitude wrt the mean Lunar surface
output.altitude = norm(yvm(1:3)) -  rm;

% Distance from the center of the Moon
output.rvm = norm(yvm(1:3));

%State vector in selenocentric coordinates
ysc = synodical2selenocentric(yv, cr3bp, tv);

%--------------------------------------------------------------------------
%Osculating circular LLO
%--------------------------------------------------------------------------
r = ysc(1:3); %position at the best approach in selenocentric frame
v = ysc(4:6); %velocity at the best approach in selenocentric frame
h =  cross(r, v); %Angular momentum vector

% Keplerian elements
[llc.a, llc.e, llc.I, llc.omega, llc.Omega] = cart2circkep(ysc);

% Mean motion
llc.n = sqrt(cr3bp.mu/llc.a^3);

% Vector of position (mean anomaly) on the LLC
Mv = 0:0.005:2*pi;
% Compute the corresponding orbit
yc = zeros(6, size(Mv, 2));
ysyn = zeros(6, size(Mv, 2));
for i = 1:size(Mv, 2)
    % Keplerian elements to selenocentric state
    yc(:,i)    = kep2cart(llc.a, llc.e, llc.I, llc.omega, llc.Omega, Mv(i), cr3bp.mu);
    % Back to synodical coordinates
    ysyn(:,i)   = selenocentric2synodical(yc(:,i), cr3bp, Mv(i)/llc.n+tv);
end

%Maneuver to perform to insert into the first orbit
output.dv1 = norm(ysc(4:6)' - yc(4:6, 1));

%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
figure(4);
hold on
quiver3(yv(1)*Lf, yv(2)*Lf, yv(3)*Lf, yv(4)*1e3*Tf, yv(5)*1e3*Tf, yv(6)*1e3*Tf, 'Color', 'b');
plot3(ysyn(1,:)*Lf, ysyn(2,:)*Lf, ysyn(3,:)*Lf, 'b');

%% Elliptic part

%We go on the llc on half an orbit
llc.T   = 2*pi/llc.n;
%Time after half an orbit:
ell.tv = llc.T/2+tv;

% Altitude of the targeted orbit: 100 km
llo.a     = 100/Lf+rm;
ell.a     = 0.5*(llc.a + llo.a);
ell.e     = (llc.a - llo.a)/(llc.a + llo.a);
ell.I     = llc.I;
ell.omega = llc.omega;
ell.Omega = llc.Omega;

% Mean motion
ell.n = sqrt(cr3bp.mu/ell.a^3);

% Vector of position (mean anomaly) on the ellptic arc
Mv = pi:0.005:2*pi;
% Compute the corresponding orbit
ell.yc = zeros(6, size(Mv, 2));
ell.ysyn = zeros(6, size(Mv, 2));
for i = 1:size(Mv, 2)
    % Keplerian elements to selenocentric state
    ell.yc(:,i)   = kep2cart(ell.a, ell.e, ell.I, ell.omega, ell.Omega, Mv(i), cr3bp.mu);
    % Back to synodical coordinates
    ell.ysyn(:,i) = selenocentric2synodical(ell.yc(:,i), cr3bp, (Mv(i)-pi)/ell.n+ell.tv);
end

%Time after half the llc + the Hohmann transfer:
llo.tv = ell.tv + 2*pi/ell.n;

%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
figure(4);
hold on
plot3(ell.ysyn(1,1)*Lf, ell.ysyn(2,1)*Lf, ell.ysyn(3,1)*Lf, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', rgb('dark green'), 'MarkerFaceColor', rgb('dark green'));
plot3(ell.ysyn(1,end)*Lf, ell.ysyn(2,end)*Lf, ell.ysyn(3,end)*Lf, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', rgb('dark blue'), 'MarkerFaceColor', rgb('dark blue'));
plot3(ell.ysyn(1,:)*Lf, ell.ysyn(2,:)*Lf, ell.ysyn(3,:)*Lf, 'b');


% Maneuver to perform to insert into the elliptic part
output.dv2 = hohmannDV1(llc.a, llo.a, cr3bp.mu);

% Maneuver to perform to insert into the llo at 50 km
output.dv3 = hohmannDV2(llc.a, llo.a, cr3bp.mu);

% Total maneuver cost
output.dvtot = output.dv1 + output.dv2 + output.dv3;

%% Targeted orbit 

llo.e = 0.0;
llo.I = llc.I;
llo.omega = llc.omega;
llo.Omega = llc.Omega;

% Mean motion
llo.n = sqrt(cr3bp.mu/llo.a^3);

% Vector of position (mean anomaly) on the LLO
Mv = 0:0.005:2*pi;
% Compute the corresponding orbit
llo.yc = zeros(6, size(Mv, 2));
llo.ysyn = zeros(6, size(Mv, 2));
for i = 1:size(Mv, 2)
    % Keplerian elements to selenocentric state
    llo.yc(:,i)    = circkep2cart(llo.a, llo.e, llo.I, llo.omega, llo.Omega, Mv(i), cr3bp.mu);
    % Back to synodical coordinates
    llo.ysyn(:,i)  = selenocentric2synodical(llo.yc(:,i), cr3bp, Mv(i)/llo.n+llo.tv);
end

%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
figure(4);
hold on
plot3(llo.ysyn(1,:)*Lf, llo.ysyn(2,:)*Lf, llo.ysyn(3,:)*Lf, 'b');
%title
strtitle = 'Reaching LLO with a Hohmann transfer';
str1     = ['DV at first injection = ', num2str(output.dv1*Tf, '%5.3f'), ' km/s'];
str2     = ['Hohmann DV1 = ', num2str(output.dv2*Tf, '%5.3f'), ' km/s'];
str3     = ['Hohmann DV2 = ', num2str(output.dv3*Tf, '%5.3f'), ' km/s'];
str4     = ['Total DV = ', num2str(output.dvtot*Tf, '%5.3f'), ' km/s'];
title({strtitle, str1, str2, str3, str4});