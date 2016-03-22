%--------------------------------------------------------------------------
% Plot a halo orbit according to settings in params.
%--------------------------------------------------------------------------
function [] = orbit_plot(yv, orbit, params, cst)
if(params.plot.XY == cst.TRUE)
    halo_orbit_plot_view(yv, 1, 2, orbit, 1);
end

if(params.plot.XZ == cst.TRUE)
    halo_orbit_plot_view(yv, 1, 3, orbit, 2);
end

if(params.plot.YZ == cst.TRUE)
    halo_orbit_plot_view(yv, 2, 3, orbit, 3);
end

if(params.plot.TD == cst.TRUE)
    halo_orbit_plot_3D(yv, orbit, 4);
end
end

%--------------------------------------------------------------------------
% Plot a halo orbit on a given plane (XY, YZ or XZ). Smallest primary
% included.
%--------------------------------------------------------------------------
function [] = halo_orbit_plot_view(yv, ip, jp, orbit, index)
%----------
%Cr3bp
%----------
cr3bp = orbit.cr3bp;

%----------
%Second primary
%----------
Rm2 = cr3bp.m2.Req/cr3bp.L;
VTheta = 0:0.01:2*pi;
%Sphere
S_L = zeros(2, size(VTheta,2));
S_L(1,:) = Rm2 *cos(VTheta);
S_L(2,:) = Rm2 *sin(VTheta);
%Position
V_L = [1 - cr3bp.mu  0  0];

%----------
%Libration point
%----------
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
plot(yv(:,ip)*Lf,yv(:,jp)*Lf, 'Color',  rgb('dark blue'), 'LineWidth', 1.5, 'LineSmoothing','on');
%Second primary
fill((V_L(ip) + S_L(1,:))*Lf,(V_L(jp) + S_L(2,:))*Lf, 'k');
%Libration point
plot(Li(ip)*Lf, Li(jp)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));

%Settings
xlabel(strcat(si, ' (x 10^3 km)'));
ylabel(strcat(sj, ' (x 10^3 km)'));
axis equal
grid on

title(strcat('Orbit projection'));
end

%--------------------------------------------------------------------------
% Plot a halo orbit in 3D. Smallest primary included.
%--------------------------------------------------------------------------
function [] = halo_orbit_plot_3D(yv, orbit, index)
%----------
%Cr3bp
%----------
cr3bp = orbit.cr3bp;

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
axis equal
grid on
%Second primary
surf(X_M3D*Lf, Y_M3D*Lf, Z_M3D*Lf, 'FaceColor', [23 153 179]./255, 'FaceLighting', 'none', 'EdgeColor', [9 63 87]./255);
%Manifold branch
plot3(yv(:,1)*Lf, yv(:,2)*Lf,yv(:,3)*Lf, 'b', 'LineWidth', 1.5, 'LineSmoothing','on');
%Libration point
plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
%Set the camera
%view([-47 28])
%zlimit
%zlim([-100 100]);
%Settings
xlabel('X (x 10^3 km)')
ylabel('Y (x 10^3 km)')
zlabel('Z (x 10^3 km)')

title(strcat('3D orbit'));
end
