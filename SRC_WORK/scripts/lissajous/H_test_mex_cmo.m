%--------------------------------------------------------------------------
% Example: this matlab file makes use of the richardson third order
% approximation to build:
% - an EML1 Lissajous orbit
% - an EML2 Lissajous orbit
%
% The Lissajous computation makes use of the mex file cmo.cpp
%
% Author: BLB
% Year: 2016
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Additionnal parameters for the cmo routine
% Order of the expansions
order = 30;
% Time span (integration)
tspan = [0 20*pi];
% Type of coordinate system as in the output
outputType = cst.coord.VSYS; 
% Type of model
model = cst.model.RTBP; %CRTBP
% Type of framework
fwrk = cst.fwrk.EM; %Earth-Moon

%% Orbit at EML1

%--------------------------------------------------------------------------
% For EML1: st0_max = 0.25 (order 40)
% For EML2: st0_max = 0.17 (order 40) 
%--------------------------------------------------------------------------
% Libration point
li = cr3bp.l1;

%--------------------------------------------------------------------------
% Initial conditions
%--------------------------------------------------------------------------
st0_max = 0.17;
st0 = st0_max*ones(1,4);

%--------------------------------------------------------------------------
%Computation with cmo
%--------------------------------------------------------------------------
[~, yo] = cmo(st0, tspan, order, li.number, outputType, model, fwrk);


%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
initplot3D(4, cr3bp, default, li)
Lf = cr3bp.L;
figure(4);
hold on
axis equal
grid on
%Manifold branch
plot3(yo(:,1)*Lf, yo(:,2)*Lf, yo(:,3)*Lf, 'Color', rgb('dark green'), 'LineWidth', 1);



%% Orbit at EML2

%--------------------------------------------------------------------------
% For EML1: st0_max = 0.25 (order 40)
% For EML2: st0_max = 0.17 (order 40) 
%-------------------------------------------------------------------------- 
% Libration point
li = cr3bp.l2;

%--------------------------------------------------------------------------
% Initial conditions
%--------------------------------------------------------------------------
st0_max = 0.12;
st0 = st0_max*ones(1,4);

%--------------------------------------------------------------------------
%Computation with cmo
%--------------------------------------------------------------------------
[~, yo] = cmo(st0, tspan, order, li.number, outputType, model, fwrk);

%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
initplot3D(4, cr3bp, default, li)
Lf = cr3bp.L;
figure(4);
hold on
axis equal
grid on
%Manifold branch
plot3(yo(:,1)*Lf, yo(:,2)*Lf, yo(:,3)*Lf, 'Color', rgb('dark blue'), 'LineWidth', 1);

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
    figure(4);
    view([-47 28]);
end