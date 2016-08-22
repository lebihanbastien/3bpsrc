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
%Constants
%--------------------------------------------------------------------------
%Distance in  km
Lf = cr3bp.L;

%--------------------------------------------------------------------------
%Create handle
%--------------------------------------------------------------------------
figure(index);
subplot(m,n,p)
hold on

%--------------------------------------------------------------------------
%Settings
%--------------------------------------------------------------------------
axis equal
grid on
zoom off
xlabel('X (x km)')
ylabel('Y (x km)')
zlabel('Z (x km)')
title ('Low Lunar Orbits (Selenocentric frame)');

%Default size
set(figure(index), 'defaultTextFontSize', 12);
set(figure(index), 'defaultTextFontWeight', 'bold');
set(figure(index), 'defaultTextHorizontalAlignment', 'center');
set(figure(index), 'defaultLineMarkerSize', 2);

%--------------------------------------------------------------------------
%Build the model
%--------------------------------------------------------------------------
moon_system(cr3bp, params);

end


%--------------------------------------------------------------------------
% Plot the Earth_Moon system
%--------------------------------------------------------------------------
function [] = moon_system(cr3bp, params)

%Distance in 10^3 km
Lf = cr3bp.L;

%----------
%Second primary
%----------
Rm2 = params.plot.bigPrimFac*cr3bp.m2.Req/cr3bp.L;
[X_M3D, Y_M3D, Z_M3D] = sphere(60);
X_M3D = 0           + Rm2*X_M3D;
Y_M3D = 0           + Rm2*Y_M3D;
Z_M3D = 0           + Rm2*Z_M3D;
HMOON = surf(X_M3D*Lf, Y_M3D*Lf, Z_M3D*Lf);

% Load Moon Image
load moonalb
% Set it on MOON
set(HMOON,'facecolor','texture',...
    'cdata',im2double(moonalb),...
    'edgecolor','none');
colormap(gray(256));

%----------
%Names
%----------
if(params.plot.names)
    if(params.plot.firstPrimDisp)
        text(-cr3bp.mu*cr3bp.L, 0,  -50, 'Earth');
    end
    text((1-cr3bp.mu)*cr3bp.L, 0, -50, 'Moon');
end


%----------
%Arrow axes
%----------
if(params.plot.tdAxes)
    
    arrow3([0 0 0], [1.3*Lf 0 0], 'k', 2, 10,0.5);
    arrow3([0 0 0], [0 0.5*Lf 0], 'k', 2, 10,0.5);
    arrow3([0 0 0], [0 0 0.5*Lf], 'k', 2, 10,0.5);
    
end

%----------
% Selenographic frame
%----------
% for i = 1:3
%   axsg = zeros(6,1);
%   axsg(i) = 1e-2;
%   c = zeros(1,3);
%   c(i) = 1;
%   arrow3(cr3bp.m2.pos*Lf, axsg(1:3)'*Lf, c, 2, 10,0.5);
% end

% Change view
view([-47 28]);

end
