function [] = initplot2D(vari, cr3bp, params, li, ip, jp)
% INITPLOT2D initializes the 2D plots used throughout the computations.
%
% INITPLOT2D(VARI, CR3BP, PARAMS, LI, IP, JP) initializes a 2D plot
% associated to the structures CR3BP (system), PARAMS (user-defined
% parameters), and LI (Lagrange point of reference). INDEX gives the number
% of the figure. The scalars IP and JP gives the dimensions that will be
% plotted. Examples:
%           - IP = 1 and JP = 2 will initialize an XY plot.
%           - IP = 2 and JP = 3 will initialize an YZ plot.
%           - IP = 3 and JP = 1 will initialize an ZX plot.
%
% If VARI is an integer, the figure initialized is figure(VARI). Else, it
% is considered that VARI is a handle to either a FIGURE or a SUBPLOT (AXES).
%
% See also INITPLOT3D
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
%String (x-axis)
switch(ip)
    case 1
        si = 'X';
    case 2
        si = 'Y';
    case 3
        si = 'Z';
    otherwise
        error('Wrong dimension on the x-axis');
end
%String (y-axis)
switch(jp)
    case 1
        sj = 'X';
    case 2
        sj = 'Y';
    case 3
        sj = 'Z';
    otherwise
        error('Wrong dimension on the y-axis');
end

xlabel(strcat(si, ' (km)'));
ylabel(strcat(sj, ' (km)'));

%--------------------------------------------------------------------------
%Constants
%--------------------------------------------------------------------------
%Distance in  10^3 km
Lf = cr3bp.L;

%--------------------------------------------------------------------------
%Second primary
%--------------------------------------------------------------------------
if(params.plot.secondPrimDisp)
    Rm2 = cr3bp.m2.Req/cr3bp.L;
    VTheta = 0:0.01:2*pi;
    %Sphere
    S_L = zeros(2, size(VTheta,2));
    S_L(1,:) = Rm2 *cos(VTheta);
    S_L(2,:) = Rm2 *sin(VTheta);
    %Position
    V_L = [1 - cr3bp.mu  0  0];
    fill((V_L(ip) + S_L(1,:))*Lf,(V_L(jp) + S_L(2,:))*Lf, 'k');
end

%--------------------------------------------------------------------------
%First primary
%--------------------------------------------------------------------------
if(params.plot.firstPrimDisp)
    
    Rm1 = cr3bp.m1.Req/cr3bp.L;
    %Sphere
    S_E = zeros(2, size(VTheta,2));
    S_E(1,:) = Rm1 *cos(VTheta);
    S_E(2,:) = Rm1 *sin(VTheta);
    %Position
    V_E = [-cr3bp.mu  0  0];
    
    fill((V_E(ip) + S_E(1,:))*Lf,(V_E(jp) + S_E(2,:))*Lf, 'k');
end

%--------------------------------------------------------------------------
%Libration point
%--------------------------------------------------------------------------
Li = li.position;
plot(Li(ip)*Lf, Li(jp)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', 2, 'MarkerSize', 2);

end
