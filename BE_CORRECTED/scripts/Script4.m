% Fourth script of the BE.
% Computing a Halo orbit around EML2 + its unstable directions + the
% complete unstable manifold, stopping at the closest approach of the lunar
% surface.
%
% Objectives: 
%  - Produce the tube-like shape of the unstable manifold.
%  - See the influence of the size of the orbit (Az) on the reachable
% altitudes with respect to the Lunar surface. Az must be taken below
% 20000 km
%
% 2016

%% Call of the first script: building the orbit
Script1

%% Stop plotting on figure 1
default.plot.XY = false;

%% Integration duration of the manifold arbitrarily fixed to 22 days
t = 22*cst.env.days*2*pi/cr3bp.T;

%% Initialization of the manifold
msi = init_manifold_branch(cst.manifold.UNSTABLE, cst.manifold.INTERIOR);

%% Initialize some useful constants
% Lunar radius
rm = cr3bp.m2.Rm/cr3bp.L;
% Lunar position
pm = cr3bp.m2.pos;
% Earth-Moon distance
Lf = cr3bp.L;

%% Manifold branch
%--------------------------------------------------------------------------
% Computation of the manifold branches that stops at the
% closest approach of the lunar surface
%--------------------------------------------------------------------------
% Vector of positions on the orbit
thetav = 0:0.01:1;
thetan = size(thetav, 2);
% Loop on all the positions on the halo orbit
li = 1;
for theta = thetav
    % Computation of the manifold branch.
    msi = manifold_branch_computation_moon(cr3bp, orbit, msi, theta, t, default, cst);
    % Save the output
    output.lfb.yv(li, :)  = msi.yve(end,1:6);
    output.lfb.tv(li)     = msi.te(end);
    output.halo.theta(li) = theta;
    li = li +1;
end

%% Postprocess

%--------------------------------------------------------------------------
% Altitude with respect to the mean Lunar surface
%--------------------------------------------------------------------------
output.altitude    = zeros(thetan, 1);
for li = 1:thetan
    %Position wrt to the Moon
    yvm = output.lfb.yv(li, 1:3) - pm;
    output.altitude(li) = norm(yvm) -  rm;
end

%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------
figure;
hold on
grid on
title('Minimum lunar altitude');
xlabel('Departure position on the halo orbit (-)');
ylabel('Minimum lunar altitude (km)');
plot(output.halo.theta, output.altitude*Lf, 'Color', 'k', 'LineWidth', 2);