%--------------------------------------------------------------------------
% Plot a manifold according to settings in params.
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------
function [] = manifold_plot(yv, orbit, manifold_branch, params, cst)

%--------------
%Color
%--------------
if(manifold_branch.stability == cst.manifold.STABLE)
    if(manifold_branch.way == cst.manifold.EXTERIOR)
        color = rgb('super green');%[51 102 0]/255;
    else
        color = rgb('super green');
    end
    
else
    if(manifold_branch.way == cst.manifold.EXTERIOR)
        color = rgb('light red');%[153 0 0]/255;
    else
        color = rgb('light red');
    end
end


if(params.plot.XY == cst.TRUE)
    manifold_plot_view(yv, 1, 2, orbit, color, 1);
end

if(params.plot.XZ == cst.TRUE)
    manifold_plot_view(yv, 1, 3, orbit,  color, 2);
end

if(params.plot.YZ == cst.TRUE)
    manifold_plot_view(yv, 2, 3, orbit, color, 3);
end

if(params.plot.TD == cst.TRUE)
    manifold_plot_3D(yv, orbit, color, 4);
end

end

%--------------------------------------------------------------------------
% Plot a manifold on a given plane (XY, YZ or XZ). Smallest primary
% included.
%--------------------------------------------------------------------------
function [] = manifold_plot_view(yv, ip, jp, orbit,  color, index)
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
    initplot2D(index, cr3bp, orbit.li, ip, jp);
end

%----------
%Plot
%----------
figure(index);
hold on
%Orbit
plot(yv(:,ip)*Lf,yv(:,jp)*Lf, 'Color',  color, 'LineWidth', 1);
end

%--------------------------------------------------------------------------
% Plot a manifold in 3D. Primaries included.
%--------------------------------------------------------------------------
function [] = manifold_plot_3D(yv, orbit,  color, index)
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

%----------
%Manifold branch
%----------
figure(index);
hold on
grid on
%Manifold branch
plot3(yv(:,1)*Lf, yv(:,2)*Lf,yv(:,3)*Lf, 'Color', color, 'LineWidth', 1);
    
end