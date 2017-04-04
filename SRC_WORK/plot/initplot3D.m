function [] = initplot3D(vari, cr3bp, params, varargin)
% INITPLOT3D initializes the 2D plots used throughout the computations.
%
% INITPLOT3D(VARI, CR3BP, PARAMS) initializes a 3D plot
% associated to the structures CR3BP (system), PARAMS (user-defined
% parameters).
%
% If VARI is an integer, the figure initialized is figure(VARI). Else, it
% is considered that VARI is a handle to either a FIGURE or a SUBPLOT (AXES).
%
% INITPLOT3D(VARI, CR3BP, PARAMS, LI) initializes a 3D plot
% associated to the structures CR3BP (system), PARAMS (user-defined
% parameters). LI must be a structure that contains one libration point.
% Kept for consistency with previous versions.
%
% INITPLOT3D(VARI, CR3BP, PARAMS, COORD) initializes a 3D plot
% associated to the structures CR3BP (system), PARAMS (user-defined
% parameters). COORD is a char that must be taken in 'syn', 'eci', 'mci',
% 'bci'.
%
% See also INITPLOT2D
%
% BLB 2016, 2017

%--------------------------------------------------------------------------
% Set figure active
%--------------------------------------------------------------------------
setplotactive(vari);
hold on;

%--------------------------------------------------------------------------
%Settings
%--------------------------------------------------------------------------
axis equal
grid on
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')

%--------------------------------------------------------------------------
%Switch between models
%--------------------------------------------------------------------------
if(size(varargin, 1) == 0)
    switch(cr3bp.name)
        case('EARTH+MOON')
            earth_moon_system(cr3bp, params);
        case('SUN+EARTH')
            sun_earth_system(cr3bp, params);
        otherwise %default is EARTH+MOON
            earth_moon_system(cr3bp, params);
    end
else
    %----------------------------------------------------------------------
    % If varargin{1} is a structure, it must be a libration point
    %----------------------------------------------------------------------
    if(isa(varargin{1}, 'struct'))
        switch(cr3bp.name)
            case('EARTH+MOON')
                earth_moon_system(cr3bp, params, varargin{1});
            case('SUN+EARTH')
                sun_earth_system(cr3bp, params,  varargin{1});
            otherwise %default is EARTH+MOON
                earth_moon_system(cr3bp, params, varargin{1});
        end
        %------------------------------------------------------------------
        % Else, we suppose that it is a char
        %------------------------------------------------------------------
    else
        switch(cr3bp.name)
            case('EARTH+MOON')
                earth_moon_system(cr3bp, params, varargin);
            case('SUN+EARTH')
                sun_earth_system(cr3bp, params); %only synodical for now
            otherwise %default is EARTH+MOON
                earth_moon_system(cr3bp, params, varargin);
        end
    end
end

end


%--------------------------------------------------------------------------
% Plot the Earth_Moon system
%--------------------------------------------------------------------------
function [] = earth_moon_system(cr3bp, params, varargin)


%--------------------------------------------------------------------------
% If there is no additional parameter, we compute the synodical frame,
% without additionnal parameter
%--------------------------------------------------------------------------
if(size(varargin, 1) == 0)
    earth_moon_syn(cr3bp, params);
elseif(isa(varargin{1}, 'struct'))
    earth_moon_syn(cr3bp, params, varargin);
else
    cst = char(varargin{1}{1});
    switch(cst)
        case('syn')
            earth_moon_syn(cr3bp, params);
        case('eci')
            earth_moon_eci(cr3bp, params);
        case('mci')
            earth_moon_mci(cr3bp, params);
        case('bci')
            earth_moon_bci(cr3bp, params);
        otherwise
            error('Unknown coordinate system');
    end
end

end

%--------------------------------------------------------------------------
% Plot the Earth_Moon system in eci coordinates
%--------------------------------------------------------------------------
function [] = earth_moon_eci(cr3bp, params)

%--------------------------------------------------------------------------
% First primary
%--------------------------------------------------------------------------
if(params.plot.firstPrimDisp)
    plot_earth(gcf, cr3bp, params, 'eci');
end

%--------------------------------------------------------------------------
% Second primary
%--------------------------------------------------------------------------
if(params.plot.secondPrimDisp)
    
    %----------------------------------------------------------------------
    % Moon position
    %----------------------------------------------------------------------
    plot_moon(cr3bp, params, 'eci');
    
    %----------------------------------------------------------------------
    % Moon's orbit
    %----------------------------------------------------------------------
    plot_primary_orbit(cr3bp, cr3bp.m2, params, 'eci');
end

end

%--------------------------------------------------------------------------
% Plot the Earth_Moon system in synodical coordinates
%--------------------------------------------------------------------------
function [] = earth_moon_syn(cr3bp, params, varargin)

%--------------------------------------------------------------------------
% Inner variables
%--------------------------------------------------------------------------
% Inner time
t0 = params.plot.tphase;

%Distance in 10^3 km
Lf = cr3bp.L;

%--------------------------------------------------------------------------
% First primary
%--------------------------------------------------------------------------
if(params.plot.firstPrimDisp)
    plot_earth(gcf, cr3bp, params, 'syn');
end

%--------------------------------------------------------------------------
% Second primary
%--------------------------------------------------------------------------
if(params.plot.secondPrimDisp)
    plot_moon(cr3bp, params, 'syn');
end

%--------------------------------------------------------------------------
%Names
%--------------------------------------------------------------------------
if(params.plot.names)
    if(params.plot.firstPrimDisp)
        text(-cr3bp.mu*cr3bp.L, 0,  -50*1e3, 'Earth');
    end
    if(params.plot.secondPrimDisp)
        text((1-cr3bp.mu)*cr3bp.L, 0, -50*1e3, 'Moon');
    end
end

%--------------------------------------------------------------------------
%Libration points, with their names on top if desired
%--------------------------------------------------------------------------
if(params.plot.allLibPoints)
    Li = cr3bp.l1.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'),  'MarkerFaceColor',  rgb('dark red'), 'MarkerSize', 5);
    if(params.plot.names)
        text(Li(1)*cr3bp.L, 0, 0.1301*Lf, 'L_1');
    end
    
    Li = cr3bp.l2.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'MarkerSize', 5);
    if(params.plot.names)
        text(Li(1)*cr3bp.L, 0, 0.1301*Lf, 'L_2');
    end
elseif(size(varargin, 1) > 0)
    li = varargin{1}{1};
    Li = li.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'),  'MarkerFaceColor',  rgb('dark red'), 'MarkerSize', 5);
    if(params.plot.names)
        text(Li(1)*cr3bp.L, 0, 0.1301*Lf, ['L_', num2str(li.number)]);
    end
end


%--------------------------------------------------------------------------
% Reference frames @t
%--------------------------------------------------------------------------
if(params.plot.frames)
    %----------------------------------------------------------------------
    % MCMF
    %----------------------------------------------------------------------
    for ii = 1:3
        c = zeros(1,3);
        c(ii) = 1;
        axsg = zeros(6,1);
        axsg(ii) = 2e-2;
        
        axsynmc =  mcmf2syn(axsg, cr3bp);
        arrow3(cr3bp.m2.pos*Lf, axsynmc(1:3)'*Lf, c, 2, 10,0.5);
    end
    
    %----------------------------------------------------------------------
    % ECEF
    %----------------------------------------------------------------------
    for ii = 1:3
        c = zeros(1,3);
        c(ii) = 1;
        axsg = zeros(6,1);
        axsg(ii) = 5e-2;
        
        axsynmc =  ecef2syn(axsg, t0, cr3bp);
        arrow3(cr3bp.m1.pos*Lf, axsynmc(1:3)'*Lf, c, 2, 10,0.5);
    end
    
    %----------------------------------------------------------------------
    % ECI
    %----------------------------------------------------------------------
    for ii = 1:3
        c = zeros(1,3);
        axsg = zeros(6,1);
        axsg(ii) = 4e-2;
        
        axsynmc =  eci2syn(axsg, t0, cr3bp);
        arrow3(cr3bp.m1.pos*Lf, axsynmc(1:3)'*Lf, c, 2, 10,0.5);
    end
    
    %----------------------------------------------------------------------
    % MCI
    %----------------------------------------------------------------------
    for ii = 1:3
        c = zeros(1,3);
        axsg = zeros(6,1);
        axsg(ii) = 1e-2;
        
        axsynmc =  mci2syn(axsg, t0, cr3bp);
        arrow3(cr3bp.m2.pos*Lf, axsynmc(1:3)'*Lf, c, 2, 10,0.5);
    end
    
    %----------------------------------------------------------------------
    % BCI, centered at the Earth, for convenience
    %----------------------------------------------------------------------
    for ii = 1:3
        c = 0.5*ones(1,3);
        axsg = zeros(6,1);
        axsg(ii) = 3e-2;
        
        axsynmc =  bci2syn(axsg, t0);
        arrow3(cr3bp.m1.pos*Lf, axsynmc(1:3)'*Lf+cr3bp.m1.pos*Lf, c, 2, 10,0.5);
    end
    
end

end

%--------------------------------------------------------------------------
% Plot the Earth_Moon system in mci coordinates
%--------------------------------------------------------------------------
function [] = earth_moon_mci(cr3bp, params)

%--------------------------------------------------------------------------
% First primary
%--------------------------------------------------------------------------
if(params.plot.firstPrimDisp)
    plot_earth(gcf, cr3bp, params, 'mci');
    
    %----------------------------------------------------------------------
    % Earth's orbit
    %----------------------------------------------------------------------
    plot_primary_orbit(cr3bp, cr3bp.m1, params, 'mci');
    
end

%--------------------------------------------------------------------------
% Second primary
%--------------------------------------------------------------------------
if(params.plot.secondPrimDisp)
    
    %----------------------------------------------------------------------
    % Moon position
    %----------------------------------------------------------------------
    plot_moon(cr3bp, params, 'mci');
end

end

%--------------------------------------------------------------------------
% Plot the Earth_Moon system in bci coordinates
%--------------------------------------------------------------------------
function [] = earth_moon_bci(cr3bp, params)

%--------------------------------------------------------------------------
% First primary
%--------------------------------------------------------------------------
if(params.plot.firstPrimDisp)
    plot_earth(gcf, cr3bp, params, 'bci');
    
    %----------------------------------------------------------------------
    % Earth's orbit
    %----------------------------------------------------------------------
    plot_primary_orbit(cr3bp, cr3bp.m1, params, 'bci');
    
end

%--------------------------------------------------------------------------
% Second primary
%--------------------------------------------------------------------------
if(params.plot.secondPrimDisp)
    
    %----------------------------------------------------------------------
    % Moon position
    %----------------------------------------------------------------------
    plot_moon(cr3bp, params, 'bci');
    
    %----------------------------------------------------------------------
    % Moon's orbit
    %----------------------------------------------------------------------
    plot_primary_orbit(cr3bp, cr3bp.m2, params, 'bci');
end

end


%--------------------------------------------------------------------------
% Plot the Moon in coordout
%--------------------------------------------------------------------------
function [] = plot_primary_orbit(cr3bp, mi, params, coordout)
    % True anomaly = Mean anomaly
    num = -pi:0.01:pi;
    % Corresponding time (the mean angular motion of the Moon is one)
    tvm = num + params.plot.tphase;
    % Orbit in synodical coordinates, then in coordout coordinates
    yvm      = zeros(6, size(num,2));
    yvm(1,:) = mi.pos(1)*ones(size(yvm(1,:)));
    yvmeci   = emcocvec(yvm, tvm, cr3bp, 'syn', coordout);
    plot3(yvmeci(1,:)*cr3bp.L, yvmeci(2,:)*cr3bp.L, yvmeci(3,:)*cr3bp.L, 'k--','LineWidth', 1);
end


%--------------------------------------------------------------------------
% Plot the Moon in coordout
%--------------------------------------------------------------------------
function [] = plot_moon(cr3bp, params, coordout)

%--------------------------------------------------------------------------
% Inner variables
%--------------------------------------------------------------------------
% Inner time
t0 = params.plot.tphase;

%Distance in 10^3 km
Lf = cr3bp.L;


%--------------------------------------------------------------------------
%Lunar globe
%--------------------------------------------------------------------------
Rm2 = params.plot.bigPrimFac*cr3bp.m2.Req/Lf;

% Sphere, in ECEF frame
[X_3D, Y_3D, Z_3D] = sphere;
X_3D = Rm2 * X_3D; Y_3D = Rm2 * Y_3D; Z_3D = Rm2 * Z_3D;


% Put X_E3D, Y_E3D, Z_E3D as rows of a matrix
xyz = [X_3D(:) Y_3D(:) Z_3D(:)].';

% MCMF to COORDOUT
for ii = 1:size(xyz,2)
    ymcmf = [xyz(1,ii), xyz(2,ii), xyz(3,ii) 0 0 0]';
    yout  = emcoc(ymcmf, t0, cr3bp, 'mcmf', coordout);
    xyz(1,ii) = yout(1);
    xyz(2,ii) = yout(2);
    xyz(3,ii) = yout(3);
end

% Reshape transformed rows
X_3D = reshape(xyz(1,:), size(X_3D));
Y_3D = reshape(xyz(2,:), size(X_3D));
Z_3D = reshape(xyz(3,:), size(X_3D));

% Plot
HMOON = surf(X_3D*Lf, Y_3D*Lf, Z_3D*Lf, 'FaceColor', [23 153 179]./255, 'FaceLighting', 'none', 'EdgeColor', [9 63 87]./255);

% Load Moon Image. CAREFUL: quite heavy for images!
load moonalb

%----------------------------------------------------------------------
% For unknown reasons, true image of the moon only works when
% the Earth is not displayed
%----------------------------------------------------------------------
if(~params.plot.firstPrimDisp)
    % Set it on MOON
    set(HMOON,'facecolor','texture','cdata', im2double(moonalb),...
        'edgecolor','none');
    colormap(gray(256));
end
end


%--------------------------------------------------------------------------
% Plot the Sun-Earth system
%--------------------------------------------------------------------------
function [] = sun_earth_system(cr3bp, params, varargin)

%Distance in km
Lf = cr3bp.L;

%--------------------------------------------------------------------------
%Arrow axes
%--------------------------------------------------------------------------
if(params.plot.tdAxes)
    arrow3([0.985*Lf 0 0], [1.015*Lf 0 0], 'k', 2, 5e-4*Lf, 0.5);
    %arrow3([0 0 0], [0 0.5*Lf 0], 'k', 2, 5e-4*Lf, 0.5);
    %arrow3([0 0 0], [0 0 0.5*Lf], 'k', 2, 5e-4*Lf, 0.5);
end

%--------------------------------------------------------------------------
%First primary
%--------------------------------------------------------------------------
if(params.plot.firstPrimDisp)
    Rm1 = params.plot.bigPrimFac*cr3bp.m1.Req/Lf;
    [X_E3D, Y_E3D, Z_E3D] = sphere;
    X_E3D = -cr3bp.mu  + Rm1*X_E3D;
    Y_E3D = 0          + Rm1*Y_E3D;
    Z_E3D = 0          + Rm1*Z_E3D;
    surf(X_E3D*Lf, Y_E3D*Lf, Z_E3D*Lf, 'FaceColor', rgb('gold'), 'FaceLighting', 'none', 'EdgeColor', 'none');
end

%--------------------------------------------------------------------------
%Second primary
%--------------------------------------------------------------------------
if(params.plot.secondPrimDisp)
    Rm2 = params.plot.bigPrimFac*cr3bp.m2.Req/Lf;
    [X_M3D, Y_M3D, Z_M3D] = sphere;
    X_M3D = 1-cr3bp.mu  + Rm2*X_M3D;
    Y_M3D = 0           + Rm2*Y_M3D;
    Z_M3D = 0           + Rm2*Z_M3D;
    surf(X_M3D*Lf, Y_M3D*Lf, Z_M3D*Lf, 'FaceColor', [0 0.6  0.9], 'FaceLighting', 'none', 'EdgeColor', 'none');
    
    %Earth globe
    %--------------------------------------------------------------------------
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
        text(-cr3bp.mu*Lf, 0,  -0.003*Lf, 'Sun');
    end
    if(params.plot.secondPrimDisp)
        text((1-cr3bp.mu)*Lf,  0,  -0.003*Lf, 'Earth');
    end
end

%--------------------------------------------------------------------------
% Moon's orbit (approximate)
%--------------------------------------------------------------------------
theta = 0:0.01:2*pi;
xM = (1-cr3bp.mu)*Lf + 480000*cos(theta);
yM = 480000*sin(theta);
zM = 0*cos(theta);
plot3(xM, yM, zM, 'k--', 'LineWidth', 2);

%--------------------------------------------------------------------------
%Libration points, with their names on top if desired
%--------------------------------------------------------------------------
if(params.plot.allLibPoints)
    Li = cr3bp.l1.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'),  'MarkerFaceColor',  rgb('dark red'), 'MarkerSize', 5);
    if(params.plot.names)
        text(Li(1)*Lf, 0, 0.005*Lf, 'L_1');
    end
    
    Li = cr3bp.l2.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'MarkerSize', 5);
    if(params.plot.names)
        text(Li(1)*Lf, 0, 0.005*Lf, 'L_2');
    end
    
    %     Li = cr3bp.l3.position;
    %     plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    %     if(params.plot.names)
    %         text(Li(1)*Lf, 0, 500, 'L_3');
    %     end
    %
    %     Li = cr3bp.l4.position;
    %     plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    %     if(params.plot.names)
    %         text(Li(1)*Lf, Li(2)*Lf, 500, 'L_4');
    %     end
    %
    %     Li = cr3bp.l5.position;
    %     plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'));
    %     if(params.plot.names)
    %         text(Li(1)*Lf, Li(2)*Lf, 500, 'L_5');
    %     end
    
elseif(size(varargin, 1) > 0)
    li = varargin{1}{1};
    Li = li.position;
    plot3(Li(1)*Lf, Li(2)*Lf, Li(3)*Lf, 'o', 'Color',  rgb('dark red'),  'MarkerFaceColor',  rgb('dark red'), 'MarkerSize', 5);
    if(params.plot.names)
        text(Li(1)*Lf, 0, 500, ['L_', num2str(li.number)]);
    end
end




end
