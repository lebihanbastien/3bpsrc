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

%--------------
%Primaries
%--------------
VTheta = 0:0.01:2*pi;
%First primary
%--------------
Rm1 = cr3bp.m1.Req/cr3bp.L;
%Sphere
S_E = zeros(2, size(VTheta,2));
S_E(1,:) = Rm1 *cos(VTheta);
S_E(2,:) = Rm1 *sin(VTheta);
%Position
V_E = [-cr3bp.mu  0  0];
%Second primary
%--------------
Rm2 = cr3bp.m2.Req/cr3bp.L;
%Sphere
S_L = zeros(2, size(VTheta,2));
S_L(1,:) = Rm2 *cos(VTheta);
S_L(2,:) = Rm2 *sin(VTheta);
%Position
V_L = [1 - cr3bp.mu  0  0];

%Libration point
Li = orbit.li.position;

%----------
%Strings
%----------
if(ip == 1 && jp == 2)
    si = 'X';
    sj = 'Y';
elseif(ip == 1 && jp == 3)
    si = 'X';
    sj = 'Z';
elseif(ip == 2 && jp == 3)
    si = 'Y';
    sj = 'Z';
end

%----------
%Factor
%----------
Lf = cr3bp.L*10^(-3);

%----------
%Plot
%----------
figure(index);
hold on
%Orbit
plot(yv(:,ip)*Lf,yv(:,jp)*Lf, 'Color',  color, 'LineWidth', 1);
%First primary
fill((V_E(ip) + S_E(1,:))*Lf,(V_E(jp) + S_E(2,:))*Lf, 'k');
%Second primary
fill((V_L(ip) + S_L(1,:))*Lf,(V_L(jp) + S_L(2,:))*Lf, 'k');
%Libration point
plot(Li(ip)*Lf, Li(jp)*Lf, 'o', 'Color',  color, 'MarkerFaceColor',  rgb('dark red'));
%Settings
xlabel(strcat(si, ' (x 10^3 km)'));
ylabel(strcat(sj, ' (x 10^3 km)'));
axis equal
grid on

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
%First primary
%----------
Rm1 = cr3bp.m1.Req/cr3bp.L;
[X_E3D, Y_E3D, Z_E3D] = sphere;
X_E3D = -cr3bp.mu  + Rm1*X_E3D;
Y_E3D = 0          + Rm1*Y_E3D;
Z_E3D = 0          + Rm1*Z_E3D;
%----------
%Second primary
%----------
Rm2 = cr3bp.m2.Req/cr3bp.L;
[X_M3D, Y_M3D, Z_M3D] = sphere;
X_M3D = 1-cr3bp.mu  + Rm2*X_M3D;
Y_M3D = 0           + Rm2*Y_M3D;
Z_M3D = 0           + Rm2*Z_M3D;


%----------
%Libration point
%----------
Li = orbit.li.position;

%----------
%Factor
%----------
Lf = cr3bp.L*10^(-3);

%----------
%Plot
%----------
figure(index);
hold on
grid on
 %First primary
surf(X_E3D*Lf, Y_E3D*Lf, Z_E3D*Lf, 'FaceColor', [23 153 179]./255, 'FaceLighting', 'none', 'EdgeColor', [9 63 87]./255);
%Second primary
surf(X_M3D*Lf, Y_M3D*Lf, Z_M3D*Lf, 'FaceColor', [23 153 179]./255, 'FaceLighting', 'none', 'EdgeColor', [9 63 87]./255);
%Manifold branch
plot3(yv(:,1)*Lf, yv(:,2)*Lf,yv(:,3)*Lf, 'Color', color, 'LineWidth', 1);
%plot3(yv(:,1)*Lf, yv(:,2)*Lf,yv(:,3)*Lf, 'LineWidth', 1);
%Libration point
plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', 2);
%Set the camera
%view([-47 28])
%Settings
xlabel('X (x 10^3 km)')
ylabel('Y (x 10^3 km)')
zlabel('Z (x 10^3 km)')
axis equal

    
end