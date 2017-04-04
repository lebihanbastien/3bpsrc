function [] = initsubplot3D(index, m, n, p, cr3bp, params, li)
% INITPLOT3D initializes the 2D plots used throughout the computations.
%
% INITPLOT3D(INDEX, CR3BP, PARAMS) initializes a 3D plot
% associated to the structures CR3BP (system), PARAMS (user-defined
% parameters). INDEX gives the number of the figure.
%
% By default, the Lagrange points L1 and L2 are plotted on the figure.
%
% See also INITPLOT2D
%
% BLB 2016

%--------------------------------------------------------------------------
%Create handle
%--------------------------------------------------------------------------
figure(index);
subplot(m,n,p)
hold on

%--------------------------------------------------------------------------
%Settings
%--------------------------------------------------------------------------
grid on
xlabel('X (x km)')
ylabel('Y (x km)')
zlabel('Z (x km)')
title ('Trajectories in Earth-Moon synodical frame');

%Default size
set(figure(index), 'defaultTextFontSize', 12);
set(figure(index), 'defaultTextFontWeight', 'bold');
set(figure(index), 'defaultTextHorizontalAlignment', 'center');
set(figure(index), 'defaultLineMarkerSize', 2);
%--------------------------------------------------------------------------
%Switch between models
%--------------------------------------------------------------------------
switch(cr3bp.name)
    case('EARTH+MOON')
        earth_moon_system(cr3bp, params, li);
    case('SUN+EARTH')
        sun_earth_system(cr3bp, params, li);
    otherwise %default is EARTH+MOON
        earth_moon_system(cr3bp, params, li);
end


end

%--------------------------------------------------------------------------
% Plot the Sun-Earth system
%--------------------------------------------------------------------------
function [] = sun_earth_system(cr3bp, params, li)

%Distance in km
Lf = cr3bp.L;

%----------
%First primary
%----------
if(params.plot.firstPrimDisp)
    Rm1 = params.plot.bigPrimFac*cr3bp.m1.Req/Lf;
    [X_E3D, Y_E3D, Z_E3D] = sphere;
    X_E3D = -cr3bp.mu  + Rm1*X_E3D;
    Y_E3D = 0          + Rm1*Y_E3D;
    Z_E3D = 0          + Rm1*Z_E3D;
    surf(X_E3D*Lf, Y_E3D*Lf, Z_E3D*Lf, 'FaceColor', rgb('gold'), 'FaceLighting', 'none', 'EdgeColor', 'none');
end

%----------
%Second primary
%----------
if(params.plot.secondPrimDisp)
    Rm2 = params.plot.bigPrimFac*cr3bp.m2.Req/Lf;
    [X_M3D, Y_M3D, Z_M3D] = sphere;
    X_M3D = 1-cr3bp.mu  + Rm2*X_M3D;
    Y_M3D = 0           + Rm2*Y_M3D;
    Z_M3D = 0           + Rm2*Z_M3D;
    surf(X_M3D*Lf, Y_M3D*Lf, Z_M3D*Lf, 'FaceColor', [0 0.6  0.9], 'FaceLighting', 'none', 'EdgeColor', 'none');
    
    %Earth globe
    %----------
    %Longitude shift
    h = 50;
    %Radius of the Earth
    r = (params.plot.bigPrimFac+0.01)*cr3bp.m2.Req;
    %Loading continents
    load coast
    lat = lat/180*pi;
    long = -(long-h)/180*pi;
    %Building the coast
    x = r*cos(lat).*sin(long);
    y = r*cos(lat).*cos(long);
    z = r*sin(lat);
    plot3((1-cr3bp.mu)*Lf +x,y,z, 'Color', [0 0 0], 'LineWidth', 1.5);
end

if(params.plot.names)
    if(params.plot.firstPrimDisp)
        text(-cr3bp.mu*Lf, 0,  -500, 'Sun');
    end
    if(params.plot.secondPrimDisp)
        text((1-cr3bp.mu)*Lf,  0,  -500, 'Earth');
    end
end

%----------
% Moon's orbit (approximate)
%----------
theta = 0:0.01:2*pi;
xM = (1-cr3bp.mu)*Lf + 480*cos(theta);
yM = 480*sin(theta);
zM = 0*cos(theta);
plot3(xM, yM, zM, 'k--', 'LineWidth', 1);

%----------
%Libration points, with their names on top if desired
%----------
if(params.plot.allLibPoints)
    Li = cr3bp.l1.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'),  'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*Lf, 0, 500, 'L_1');
    end
    
    Li = cr3bp.l2.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*Lf, 0, 500, 'L_2');
    end
    
    Li = cr3bp.l3.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*Lf, 0, 500, 'L_3');
    end
    
    Li = cr3bp.l4.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*Lf, Li(2)*Lf, 500, 'L_4');
    end
    
    Li = cr3bp.l5.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*Lf, Li(2)*Lf, 500, 'L_5');
    end
    
else
    Li = li.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'),  'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*Lf, 0, 500, ['L_', num2str(li.number)]);
    end
end


%----------
%Arrow axes
%----------
if(params.plot.tdAxes)
    arrow3([1.475e5 0 0], [1.52e5 0 0], 'k', 2, 50,0.5);
end

end

%--------------------------------------------------------------------------
% Plot the Earth_Moon system
%--------------------------------------------------------------------------
function [] = earth_moon_system(cr3bp, params, li)

%Distance in 10^3 km
Lf = cr3bp.L;

%----------
%First primary
%----------
if(params.plot.firstPrimDisp)
    Rm1 = params.plot.bigPrimFac*cr3bp.m1.Req/cr3bp.L;
    [X_E3D, Y_E3D, Z_E3D] = sphere;
    X_E3D = -cr3bp.mu  + Rm1*X_E3D;
    Y_E3D = 0          + Rm1*Y_E3D;
    Z_E3D = 0          + Rm1*Z_E3D;
    surf(X_E3D*Lf, Y_E3D*Lf, Z_E3D*Lf, 'FaceColor', [0 0.6  0.9], 'FaceLighting', 'none', 'EdgeColor', 'none');
    
    %Earth globe
    %----------
    %Longitude shift
    h = 50;
    %Radius of the Earth
    r = (params.plot.bigPrimFac+0.01)*cr3bp.m1.Req;
    %Loading continents
    load coast
    lat = lat/180*pi;
    long = -(long-h)/180*pi;
    %Building coast
    x = r*cos(lat).*sin(long);
    y = r*cos(lat).*cos(long);
    z = r*sin(lat);
    plot3(-cr3bp.mu*Lf+x,y,z, 'Color', [0 0 0 ], 'LineWidth', 1);
    
end

%----------
%Second primary
%----------
if(params.plot.secondPrimDisp)
    Rm2 = params.plot.bigPrimFac*cr3bp.m2.Req/cr3bp.L;
    [X_M3D, Y_M3D, Z_M3D] = sphere;
    X_M3D = 1-cr3bp.mu  + Rm2*X_M3D;
    Y_M3D = 0           + Rm2*Y_M3D;
    Z_M3D = 0           + Rm2*Z_M3D;
    surf(X_M3D*Lf, Y_M3D*Lf, Z_M3D*Lf, 'FaceColor', [23 153 179]./255, 'FaceLighting', 'none', 'EdgeColor', [9 63 87]./255);
end

%----------
%Names
%----------
if(params.plot.names)
    if(params.plot.firstPrimDisp)
        text(-cr3bp.mu*cr3bp.L, 0,  -50, 'Earth');
    end
    if(params.plot.secondPrimDisp)
        text((1-cr3bp.mu)*cr3bp.L, 0, -50, 'Moon');
    end
end

%----------
%Libration points, with their names on top if desired
%----------
if(params.plot.allLibPoints)
    Li = cr3bp.l1.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'),  'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*cr3bp.L, 0, 0.1301*Lf, 'L_1');
    end
    
    Li = cr3bp.l2.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*cr3bp.L, 0, 0.1301*Lf, 'L_2');
    end
    
    Li = cr3bp.l3.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*cr3bp.L, 0, 0.0780*Lf, 'L_3');
    end
    
    Li = cr3bp.l4.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*cr3bp.L, Li(2)*cr3bp.L, 0.0780*Lf, 'L_4');
    end
    
    
    Li = cr3bp.l5.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*cr3bp.L, Li(2)*cr3bp.L, 0.0780*Lf, 'L_5');
    end
    
else
    Li = li.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'),  'MarkerFaceColor',  rgb('dark red'));
    if(params.plot.names)
        text(Li(1)*cr3bp.L, 0, 0.1301*Lf, ['L_', num2str(li.number)]);
    end
end


%----------
%Arrow axes
%----------
if(params.plot.tdAxes)
    
    arrow3([0 0 0], [1.3*Lf 0 0], 'k', 2, 10,0.5);
    arrow3([0 0 0], [0 0.5*Lf 0], 'k', 2, 10,0.5);
    arrow3([0 0 0], [0 0 0.5*Lf], 'k', 2, 10,0.5);
    
end
end
