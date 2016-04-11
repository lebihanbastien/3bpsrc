%--------------------------------------------------------------------------
% Example n°1: this matlab file makes use of the richardson third order
% approximation to build:
% - an EML2 planar lyapunov orbit (~8000 km of radius)
% - an EML2 halo orbit (~10 000km of vertical extension)
% - an EML2 vertical lyapunov orbit (~30 000km of vertical extension)
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Environment init
cr3bp = init_CR3BP('SUN', 'EARTH', default);

%% Init plot
initplot3D(4, cr3bp, default);


%% Test of cumovec.cpp

% For CRTBP EML1: st0_max = 0.25 (order 40)
% For CRTBP EML2: st0_max = 0.17 (order 40)

% @TODO: 1. CRTBP compute up to order 40 again (currently 30)
%        2. Check the order for QBCP (currently unknown, probably 20).
%        3. Add the event option: stop integration if the solution is going
%        to far away from SEML1,2
%        4. Make the inputs more robust: no possibility to have kmax >
%        ngrid for example: segmentation fault is around the corner!
%        Moreover, externalise ALL the number of points: along the orbit +
%        along the manifold... make it flexible! Maybe a time vector
%        instead of [kmin, kmax] + kstep?
%        5. The time vector tspan should be given in terms of either orbit
%        period or Sun period. Sun period is good, because that way we can
%        isolate more easily the "good" regions. In particular, look if we
%        can isolate the good regions in terms of timespan
%        6. Think of dichotomy to isolate the "good" solutions that stay
%        around SEML1,2
%        8. Plotter l'énergie absolue autour de SEMLi et le long des
%        variétés instables partant de EMLi.
%        9. Mettre la position de la Lune en sortie sur yMan?
%        10. Réfléchir à l'ordre des variables dans les tenseurs et
%        matrices... Pas forcément optimal actuellement pour MATLAB!
%
%
%
%        A. Change the computation: if output is just tf, yf, do NOT get
%        the results on a grid! Extra work for nothing!
%        B. Do the SEML case!

%----------
% Parameters
%----------
% Initial conditions
st0_max = 6.0;
st0     = st0_max*ones(1,4);
st0(5)  = 0.0;

% Libration point
li = cr3bp.l2;
% Order of the expansions
order = 12;
% Time span (in Sun period)
tspan = [0 1];
% Time on the manifold (in Sun period)
tman = 10;

% Inputs
outputType = cst.coord.VSEM;
model      = cst.model.QBCP;
fwrk       = cst.fwrk.EM;
manType    = cst.mantype.MAN_CENTER_U;

%Steps along the orbit
kmin  = 0;
kmax  = 50;
kstep = 2;

%----------
%Computation
%----------
[to, yo, tm, ym] = cumovec(st0, tspan, tman, order, li.number, outputType, model, fwrk, manType, kmin, kmax, kstep);


%% Plot loop
Lf = cr3bp.L;

figure(4);
hold on
grid on
%Orbit branch
plot3(yo(:,1)*Lf, yo(:,2)*Lf, yo(:,3)*Lf, 'Color', rgb('dark green'), 'LineWidth', 2);

%Manifold branch
for k = 1:2:size(ym, 3)
    %Find the last non zero value
    ind  = ym(:,1,k) ~= 0;
    yvec = ym(ind, :, k); 
    %Condition on the last state
    dist2Origin = norm(yvec(end,1))*Lf;
    if(dist2Origin < 1.52e8 && dist2Origin > 1.47e8)
        plot3(yvec(:,1)*Lf, yvec(:,2)*Lf, yvec(:,3)*Lf, 'Color', rgb('dark red'), 'LineWidth', 2);
    end
end

% Change axis style
axis square
%% Plot size of manifold around SEML2

% Get data
filename = '~/BackUpBox/PhD/OOFTDA/fprint/QBCP/SEM/L2/eIm_rand_ofs_30_order_20.bin';
fileID = fopen(filename);
Mroot  = fread(fileID, [20, Inf], 'double');
fclose(fileID);
Mroot = Mroot';

% Select only some values
M = Mroot;

% Plot
figure(4);
hold on;
grid on;
% Plot the whole data set
hp = scatter3(-M(:,8)*Lf, -M(:,9)*Lf, M(:,10)*Lf, 20, M(:,20), 'filled');
%hp = scatter3(-M(:,8)*Lf, -M(:,9)*Lf, 0*M(:,10)*Lf, 20, M(:,20), 'filled');
colormap(jet(5));
cp = colorbar;
caxis([0,5e-6])

%% Plot size of manifold around SEML1

% Get data
filename = '~/BackUpBox/PhD/OOFTDA/fprint/QBCP/SEM/L1/eIm_rand_ofs_30_order_15.bin';
fileID = fopen(filename);
Mroot  = fread(fileID, [20, Inf], 'double');
fclose(fileID);
Mroot = Mroot';

% Select only some values
M = Mroot;

% Plot
figure(4);
hold on;
grid on;
% Plot the whole data set
hp = scatter3(-M(:,8)*Lf, -M(:,9)*Lf, 0*M(:,10)*Lf, 15, M(:,20), 'filled');
%hp = scatter3(-M(:,8)*Lf, -M(:,9)*Lf, 0*M(:,10)*Lf, 20, M(:,20), 'filled');
colormap(jet(10));
cp = colorbar;
ylabel(cp, 'Invariance error');
%set(cp, 'YAxisLocation', 'left');
caxis([0,1e-5])

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
    figure(4);
    view([-47 28]);
end