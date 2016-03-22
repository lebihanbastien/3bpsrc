%--------------------------------------------------------------------------
% Plot an orbit according to settings in params.
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------
function [] = orbit_plot(yv, orbit, params, cst, isOrbitOnly)
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
    if(strcmp(orbit.cr3bp.name ,'SUN+EARTH'))
        halo_orbit_plot_3D(yv, orbit, 4, isOrbitOnly, 0, 0, 1, 1, 1);
    else
        halo_orbit_plot_3D(yv, orbit, 4, isOrbitOnly, 1, 1, 1, 1, 0);
    end
end
end

%--------------------------------------------------------------------------
% Plot an orbit on a given plane (XY, YZ or XZ). Smallest primary
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
plot(yv(:,ip)*Lf,yv(:,jp)*Lf, 'Color',  rgb('dark blue'), 'LineWidth', 1.5);
%plot(yv(:,ip)*Lf,yv(:,jp)*Lf, 'LineWidth', 1.5);
%Second primary
fill((V_L(ip) + S_L(1,:))*Lf,(V_L(jp) + S_L(2,:))*Lf, 'k');
%Libration point
Li = orbit.cr3bp.l1.position;
plot(Li(ip)*Lf, Li(jp)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', 2, 'MarkerSize', 2);

Li = orbit.cr3bp.l2.position;
plot(Li(ip)*Lf, Li(jp)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', 2, 'MarkerSize', 2);

Li = orbit.cr3bp.l3.position;
plot(Li(ip)*Lf, Li(jp)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', 2, 'MarkerSize', 2);

Li = orbit.cr3bp.l4.position;
plot(Li(ip)*Lf, Li(jp)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', 2, 'MarkerSize', 2);

Li = orbit.cr3bp.l5.position;
plot(Li(ip)*Lf, Li(jp)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', 2, 'MarkerSize', 2);
%Settings
xlabel(strcat(si, ' (x 10^3 km)'));
ylabel(strcat(sj, ' (x 10^3 km)'));
axis equal
grid on

title(strcat('Orbit projection'));
end

%--------------------------------------------------------------------------
% Plot an orbit in 3D. Smallest primary included.
%--------------------------------------------------------------------------
function [] = halo_orbit_plot_3D(yv, orbit, index, isOrbitOnly, isFirstPrimary, isAllLibPoints, isNames, isAxes, isMoonOrbit)

%----------
%Cr3bp
%----------
cr3bp = orbit.cr3bp;

%----------
%First primary
%----------
Rm1 = 2*cr3bp.m1.Req/cr3bp.L;
[X_E3D, Y_E3D, Z_E3D] = sphere;
X_E3D = -cr3bp.mu  + Rm1*X_E3D;
Y_E3D = 0          + Rm1*Y_E3D;
Z_E3D = 0          + Rm1*Z_E3D;


%----------
%Second primary
%----------
Rm2 = 5*cr3bp.m2.Req/cr3bp.L;
[X_M3D, Y_M3D, Z_M3D] = sphere;
X_M3D = 1-cr3bp.mu  + Rm2*X_M3D;
Y_M3D = 0           + Rm2*Y_M3D;
Z_M3D = 0           + Rm2*Z_M3D;

%----------
%Factor
%----------
Lf = cr3bp.L*10^(-3);

%----------
%Text
%----------
t1 = zeros(7,1);
t1index = 1;

%----------
%Plot
%----------
figure(index);
hold on

%Settings
axis equal
grid on
xlabel('X (x 10^3 km)')
ylabel('Y (x 10^3 km)')
zlabel('Z (x 10^3 km)')

%Orbit
plot3(yv(:,1)*Lf, yv(:,2)*Lf,yv(:,3)*Lf, 'Color',  rgb('dark blue'), 'LineWidth', 1.5);

if(~isOrbitOnly)
    switch(cr3bp.name)
        
        case('EARTH+MOON')
            %First primary (Earth)
            if(isFirstPrimary)
                surf(X_E3D*Lf, Y_E3D*Lf, Z_E3D*Lf, 'FaceColor', [0 0.6  0.9], 'FaceLighting', 'none', 'EdgeColor', 'none');
            end
            
            %Second primary (Moon)
            surf(X_M3D*Lf, Y_M3D*Lf, Z_M3D*Lf, 'FaceColor', [23 153 179]./255, 'FaceLighting', 'none', 'EdgeColor', [9 63 87]./255);
            
            if(isFirstPrimary)
                %Earth globe
                h = 50; %longitude shift
                r = 2.01*cr3bp.m1.Req*10^(-3); %radius of the Earth
                load coast %loading continents
                lat = lat/180*pi;
                long = -(long-h)/180*pi;
                x = r*cos(lat).*sin(long);
                y = r*cos(lat).*cos(long);
                z = r*sin(lat);
                cote = plot3(-cr3bp.mu*Lf+x,y,z);
                cote.Color = 0*[1 1 1];
                cote.LineWidth = 1;
            end
            
            if(isNames)
                if(isFirstPrimary)
                    t1(t1index) = text(-cr3bp.mu*cr3bp.L*10^(-3), 0,  -50, 'Earth');
                    t1index = t1index + 1;
                end
                t1(t1index) = text((1-cr3bp.mu)*cr3bp.L*10^(-3), 0, -50, 'Moon');
                t1index = t1index + 1;
            end
            
            
            
        case('SUN+EARTH')
            %First primary (Sun)
            if(isFirstPrimary)
                surf(X_E3D*Lf, Y_E3D*Lf, Z_E3D*Lf, 'FaceColor', rgb('gold'), 'FaceLighting', 'none', 'EdgeColor', 'none');
            end
            %Second primary (Earth)
            surf(X_M3D*Lf, Y_M3D*Lf, Z_M3D*Lf, 'FaceColor', [0 0.6  0.9], 'FaceLighting', 'none', 'EdgeColor', 'none');
            
            % Earth globe
            h = 50; %longitude shift
            r = 5.01*cr3bp.m2.Req*10^(-3); %radius of the Earth
            load coast %loading continents
            lat = lat/180*pi;
            long = -(long-h)/180*pi;
            x = r*cos(lat).*sin(long);
            y = r*cos(lat).*cos(long);
            z = r*sin(lat);
            cote = plot3((1-cr3bp.mu)*Lf +x,y,z);
            cote.Color = 0*[1 1 1];
            cote.LineWidth = 1;
            
            if(isNames)
                if(isFirstPrimary)
                    t1(t1index) = text(-cr3bp.mu*cr3bp.L*10^(-3), 0,  -5e2, 'Sun');
                    t1index = t1index + 1;
                end
                t1(t1index) = text((1-cr3bp.mu)*cr3bp.L*10^(-3), 0,  -5e2, 'Earth');
                t1index = t1index + 1;
            end
            
            
        otherwise %default is EARTH+MOON
            %First primary (Earth)
            if(isFirstPrimary)
                surf(X_E3D*Lf, Y_E3D*Lf, Z_E3D*Lf, 'FaceColor', [0 0.6  0.9], 'FaceLighting', 'none', 'EdgeColor', 'none');
            end
            %Second primary (Moon)
            surf(X_M3D*Lf, Y_M3D*Lf, Z_M3D*Lf, 'FaceColor', [23 153 179]./255, 'FaceLighting', 'none', 'EdgeColor', [9 63 87]./255);
            
            if(isFirstPrimary)
                %Earth globe
                h = 50; %longitude shift
                r = 2.01*cr3bp.m1.Req*10^(-3); %radius of the Earth
                load coast %loading continents
                lat = lat/180*pi;
                long = -(long-h)/180*pi;
                x = r*cos(lat).*sin(long);
                y = r*cos(lat).*cos(long);
                z = r*sin(lat);
                cote = plot3(-cr3bp.mu*Lf+x,y,z);
                cote.Color = 0*[1 1 1];
                cote.LineWidth = 1;
            end
            
            if(isNames)
                if(isFirstPrimary)
                    t1(t1index) = text(-cr3bp.mu*cr3bp.L*10^(-3), 0,  -50, 'Earth');
                    t1index = t1index + 1;
                end
                t1(t1index) = text((1-cr3bp.mu)*cr3bp.L*10^(-3), 0, -50, 'Moon');
                t1index = t1index + 1;
            end
    end
    
    %----------
    %Libration points, with their names on top if desired
    %----------
    psize = 5;
    if(isAllLibPoints)
        Li = orbit.cr3bp.l1.position;
        plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', psize);
        
        Li = orbit.cr3bp.l2.position;
        plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', psize);
        
        Li = orbit.cr3bp.l3.position;
        plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', psize);
        
        Li = orbit.cr3bp.l4.position;
        plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', psize);
        
        Li = orbit.cr3bp.l5.position;
        plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', psize);
        
        
        if(isNames)
            Li = orbit.cr3bp.l1.position;
            t1(t1index) = text(Li(1)*cr3bp.L*10^(-3), 0, 0.1301*Lf, 'L_1');
            t1index = t1index + 1;
            
            Li = orbit.cr3bp.l2.position;
            t1(t1index) = text(Li(1)*cr3bp.L*10^(-3), 0, 0.1301*Lf, 'L_2');
            t1index = t1index + 1;
            
            Li = orbit.cr3bp.l3.position;
            t1(t1index) = text(Li(1)*cr3bp.L*10^(-3), 0, 0.0780*Lf, 'L_3');
            t1index = t1index + 1;
            
            Li = orbit.cr3bp.l4.position;
            t1(t1index) = text(Li(1)*cr3bp.L*10^(-3), Li(2)*cr3bp.L*10^(-3), 0.0780*Lf, 'L_4');
            t1index = t1index + 1;
            
            Li = orbit.cr3bp.l5.position;
            t1(t1index) = text(Li(1)*cr3bp.L*10^(-3), Li(2)*cr3bp.L*10^(-3), 0.0780*Lf, 'L_5');
        end
        
        
    else
        Li = orbit.li.position;
        plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', psize);
        
        if(isNames)
            t1(t1index) = text(Li(1)*cr3bp.L*10^(-3), Li(2)*cr3bp.L*10^(-3), 5e2, ['L_', num2str(orbit.li.number)]);
        end
    end
    
    %Right font for the names
    if(isNames)
        for i = 1:t1index
            set(t1(i), 'FontSize', 16);
            set(t1(i), 'FontWeight', 'Bold');
            set(t1(i), 'HorizontalAlignment', 'center');
        end
    end
    
    %----------
    %Arrow axes
    %----------
    if(isAxes)
        switch(cr3bp.name)
            case('EARTH+MOON')
                arrow3([0 0 0], [1.3*Lf 0 0], 'k', 2, 10,0.5);
                arrow3([0 0 0], [0 0.5*Lf 0], 'k', 2, 10,0.5);
                arrow3([0 0 0], [0 0 0.5*Lf], 'k', 2, 10,0.5);
                
            case('SUN+EARTH')
                arrow3([1.475e5 0 0], [1.52e5 0 0], 'k', 2, 50,0.5);
        end
    end
    
    %----------
    %Moon's orbit
    %----------
    if(isMoonOrbit)
        theta = 0:0.01:2*pi;
        xM = (1-cr3bp.mu)*Lf + 480*cos(theta);
        yM = 480*sin(theta);
        zM = 0*cos(theta);
        plot3(xM, yM, zM, 'k--', 'LineWidth', 1);
    end
end
%title(strcat('3D orbit'));
end
