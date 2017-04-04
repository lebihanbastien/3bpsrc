function [ output ] = moonreachability(cr3bp, orbit, msi, thetav, maxalt, params, cst)
% OUTPUT = MOONREACHIBILITY( CR3BP, ORBIT, MSI, THETAV, MAXALT, PARAMS, CST)
% studies the possibilities of lunar satellisation from a
% selected set of positions THETAV on the orbit ORBIT and transfering on
% the manifold of type MSI. The study is limited to the manifold branches
% whose closest approach of the lunar surface is below the altitude MAXALT,
% given in normalized units.
% 

%--------------------------------------------------------------------------
% Integration duration of the manifold arbitrarily fixed to 20 days
%--------------------------------------------------------------------------
t = 20*cst.env.days*2*pi/cr3bp.T;

%--------------------------------------------------------------------------
% Computation of the manifold branches that stops at the
% closest approach of the lunar surface
%--------------------------------------------------------------------------
thetan = size(thetav, 2);

% Loop on all the positions on the halo orbit
li = 1;
for theta = thetav
    % Computation of the manifold branch.
    msi = manifold_branch_computation_moon(cr3bp, orbit, msi, theta, t, params, cst);
    % Save the output
    output.lfb.yv(li, :)  = msi.yve(end,1:6);
    output.lfb.tv(li)     = msi.te(end);
    output.halo.theta(li) = theta;
    li = li +1;
end

%--------------------------------------------------------------------------
% Initialize some useful constants
%--------------------------------------------------------------------------
% Lunar radius
rm = cr3bp.m2.Rm/cr3bp.L;
% Lunar position
pm = cr3bp.m2.pos;
% Earth-Moon distance
Lf = cr3bp.L;
% 1000 x the velocity factor
Tf = 1e3*cr3bp.L*2*pi/cr3bp.T;

%--------------------------------------------------------------------------
% Initialize the outputs
%--------------------------------------------------------------------------
% Vector of position (mean anomaly) on the LLO
Mv = 0:0.01:2*pi;
output.altitude    = zeros(thetan, 1);
output.latitude    = zeros(thetan, 1);
output.longitude   = zeros(thetan, 1);
output.llo.lat     = zeros(thetan, size(Mv, 2));
output.llo.long    = zeros(thetan, size(Mv, 2));

%--------------------------------------------------------------------------
% Latitude and Longitude on the Moon
%--------------------------------------------------------------------------
f1 = figure(6);
set(f1, 'Position', [20, 1, 1920, 959]);
subplot(2,3,1:2);
hold on
load moonalb
geoshow(moonalb, moonalbrefvec)
grid on
title('Selenographic frame');

%--------------------------------------------------------------------------
% Final Low Lunar orbits
%--------------------------------------------------------------------------
initsubplotMoon(6, 2, 3, 3, cr3bp, params, cr3bp.l2);

%--------------------------------------------------------------------------
% Loop on all the positions on the Halo orbit
%--------------------------------------------------------------------------
for li = 1:thetan
    % Current state
    yv = output.lfb.yv(li, :)';
    tv = output.lfb.tv(li);
    
    % Current state wrt the Moon
    yvm = yv - [pm 0 0 0]';
    
    % Altitude wrt the mean Lunar surface
    output.altitude(li) = norm(yvm(1:3)) -  rm;
    
    % Distance from the center of the Moon
    output.rvm(li) = norm(yvm(1:3));
    
    % Latitude and longitude in selenographic coordinates
    ysg = synodical2selenographic(yv, cr3bp);
    output.latitude(li) = rad2deg(atan2(ysg(3), norm(ysg(1:2))));
    output.longitude(li)= rad2deg(atan2(ysg(2), ysg(1)));
    
    %State vector in selenocentric coordinates
    ysc = synodical2selenocentric(yv, cr3bp, tv);
    
    
    %----------------------------------------------------------------------
    %Osculating circular LLO
    %----------------------------------------------------------------------
    r = ysc(1:3); %position at the best approach in selenocentric frame
    v = ysc(4:6); %velocity at the best approach in selenocentric frame
    h =  cross(r, v); %Angular momentum vector
    
    % Keplerian elements
    [a, e, I, omega, Omega] = cart2circkep(ysc);
    
    % Mean motion
    n = sqrt(cr3bp.mu/a^3);
    
    % Compute the corresponding orbit
    yc = zeros(6, size(Mv, 2));
    ysyn = zeros(6, size(Mv, 2));
    for i = 1:size(Mv, 2)
        % Keplerian elements to selenocentric state
        yc(:,i)    = circkep2cart(a, e, I, omega, Omega, Mv(i), cr3bp.mu);
        % Back to synodical coordinates
        ysg         = selenocentric2selenographic(yc(:,i), cr3bp, Mv(i)/n+tv);
        ysyn(:,i)   = selenocentric2synodical(yc(:,i), cr3bp, Mv(i)/n+tv);
        % Store in output
        output.llo.lat(li, i)  = rad2deg(atan2(ysg(3), norm(ysg(1:2))));
        output.llo.long(li, i) = rad2deg(atan2(ysg(2), ysg(1)));
    end
    
    % Maneuver to perform to insert into the LLO
    output.llo.dv(li) = norm(ysc(4:6)' - yc(4:6, 1));
    
    
    %----------------------------------------------------------------------
    %Plot solutions if they are good enough
    %----------------------------------------------------------------------
    if(output.altitude(li) > 0 && output.altitude(li) < maxalt)
        
        % Recompute the manifold, stopping at the right event, with plot
        params.plot.manifold_branch = true;
        msi  = manifold_branch_computation_moon(cr3bp, orbit, msi, output.halo.theta(li), t, params, cst);
        
        %Equivalent of the trajectory in selenographic coordinates
        ntraj = size(msi.ytraj, 1);
        lat = zeros(1, ntraj);
        long = zeros(1, ntraj);
        for i = 1:ntraj
            yseltraj = synodical2selenographic(msi.ytraj(i,:)', cr3bp);
            lat(i)  = rad2deg(atan2(yseltraj(3), norm(yseltraj(1:2))));
            long(i) = rad2deg(atan2(yseltraj(2), yseltraj(1)));
        end
        
        %------------------------------------------------------------------
        % Update figure 5 & 6
        %------------------------------------------------------------------
        % Synodical frame
        figure(4);
        hold on
        quiver3(yv(1)*Lf, yv(2)*Lf, yv(3)*Lf, yv(4)*Tf, yv(5)*Tf, yv(6)*Tf, 'Color', 'g');
        plot3(ysyn(1,1)*Lf, ysyn(2,1)*Lf, ysyn(3,1)*Lf, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
        plot3(ysyn(1,:)*Lf, ysyn(2,:)*Lf, ysyn(3,:)*Lf, 'Color', 'g');
        
        figure(6);
        subplot(2,3,3);
        hold on
        %quiver3(ysc(1)*Lf, ysc(2)*Lf, ysc(3)*Lf, ysc(4)*Tf, ysc(5)*Tf, ysc(6)*Tf, 'Color', 'b');
        plot3(yc(1,1)*Lf, yc(2,1)*Lf, yc(3,1)*Lf, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
        plot3(yc(1,:)*Lf, yc(2,:)*Lf, yc(3,:)*Lf, 'Color', 'g');
        zoom off
        
        % Selenographic frame
        figure(6);
        subplot(2,3,1:2);
        hold on
        geoshow(lat, long, 'DisplayType', 'Line', 'Color', 'r');
        geoshow(output.latitude(li), output.longitude(li), 'DisplayType', 'Point', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
        geoshow(output.llo.lat(li,:), output.llo.long(li,:), 'DisplayType', 'Line', 'Color', 'g');
        xlabel('Longitude (°)');
        ylabel('Latitude (°)');
        
        % Lunar altitude VS position on the initial halo orbit
        subplot(2,3,4);
        hold on
        plot(output.halo.theta(li), output.lfb.tv(li)*cr3bp.T/(2*pi*cst.env.days), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        grid on
        title('Arrival from Halo orbit');
        xlabel('Departure position on the halo orbit (-)');
        ylabel('Time of Flight (days)');
        
        % Inclination vs Altitude of the LLO
        subplot(2,3,6);
        hold on;
        xlabel('Inclination (°)');
        ylabel('Altitude (km)');
        plot(rad2deg(I), output.altitude(li)*Lf, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        grid on;
        title('Parameters of the Low Lunar Orbit (LLO)');
        
        % Inclination vs DV of the LLO
        subplot(2,3,5);
        hold on;
        plot(rad2deg(I), output.llo.dv(li)*Lf*2*pi/(cr3bp.T), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        xlabel('Inclination (°)');
        ylabel('DV (km/s)');
        grid on;
        title('Lunar insertion maneuver cost');
        
        
    end
    
end
end

