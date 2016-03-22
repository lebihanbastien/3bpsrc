%--------------------------------------------------------------------------
% Init the 2D plot (primaries...)
%--------------------------------------------------------------------------
function [] = initplot2D(index, cr3bp, params, li, ip, jp)

%--------------------------------------------------------------------------
%Constants
%--------------------------------------------------------------------------
%Distance in  10^3 km
Lf = cr3bp.L*10^(-3);

%--------------------------------------------------------------------------
%Create handle
%--------------------------------------------------------------------------
figure(index);
hold on

%--------------------------------------------------------------------------
%Settings
%--------------------------------------------------------------------------
axis equal
grid on
%Strings
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
xlabel(strcat(si, ' (x 10^3 km)'));
ylabel(strcat(sj, ' (x 10^3 km)'));
title(['Projection in the ', si, sj, '-plane']);
%Default size
set(figure(index), 'defaultTextFontSize', 12);
set(figure(index), 'defaultTextFontWeight', 'bold');
set(figure(index), 'defaultTextHorizontalAlignment', 'center');
set(figure(index), 'defaultLineMarkerSize', 2);

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
fill((V_L(ip) + S_L(1,:))*Lf,(V_L(jp) + S_L(2,:))*Lf, 'k');

%----------
%First primary
%----------
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

%----------
%Libration point
%----------
Li = li.position;
plot(Li(ip)*Lf, Li(jp)*Lf, 'o', 'Color',  rgb('dark red'), 'MarkerFaceColor',  rgb('dark red'), 'markers', 2, 'MarkerSize', 2);

%----------
%Second primary
%----------



end
