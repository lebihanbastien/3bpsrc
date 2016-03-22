%--------------------------------------------------------------------------
% Plot an orbit according to settings in params.
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------
function [] = orbit_plot(yv, orbit, params, cst)
if(params.plot.XY == cst.TRUE)
    halo_orbit_plot_view(yv, 1, 2, orbit, 1, params);
end

if(params.plot.XZ == cst.TRUE)
    halo_orbit_plot_view(yv, 1, 3, orbit, 2, params);
end

if(params.plot.YZ == cst.TRUE)
    halo_orbit_plot_view(yv, 2, 3, orbit, 3, params);
end

if(params.plot.TD == cst.TRUE)
    halo_orbit_plot_3D(yv, orbit, 4, params);
end
end

%--------------------------------------------------------------------------
% Plot an orbit on a given plane (XY, YZ or XZ). Smallest primary
% included
%--------------------------------------------------------------------------
function [] = halo_orbit_plot_view(yv, ip, jp, orbit, index, params)
%----------
%Cr3bp
%----------
cr3bp = orbit.cr3bp;

%----------
%Factor
%----------
Lf = cr3bp.L*10^(-3);

%----------
%Plot
%----------
%If the figure did not exist before, we set the basic environment.
if(~ishandle(index))
    initplot2D(index, cr3bp, params, orbit.li, ip, jp);
end

%----------
%Orbit
%----------
figure(index);
%Orbit
plot(yv(:,ip)*Lf,yv(:,jp)*Lf, 'Color',  rgb('dark blue'), 'LineWidth', 1.5);

end

%--------------------------------------------------------------------------
% Plot an orbit in 3D. Various options are available through params
%--------------------------------------------------------------------------
function [] = halo_orbit_plot_3D(yv, orbit, index, params)

%----------
%Cr3bp
%----------
cr3bp = orbit.cr3bp;

%----------
%Factor
%----------
Lf = cr3bp.L*10^(-3);

%----------
%Plot
%----------
%If the figure did not exist before, we set the basic environment.
if(~ishandle(index))
    initplot3D(index, cr3bp, params);
end

%Orbit
figure(index);
plot3(yv(:,1)*Lf, yv(:,2)*Lf,yv(:,3)*Lf, 'Color',  rgb('dark blue'), 'LineWidth', 1.5);

end

