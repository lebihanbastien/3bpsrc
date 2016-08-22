init;

%--------------------------------------------------------------------------
% Test of FMINCON for QBCP
%--------------------------------------------------------------------------
% Number of points on each legs
N1 = 100;
N2 = 3;
N3 = 100;

inputType  = cst.coord2.PSEM;
outputType = cst.coord2.VNCSEM;
dcs        = cst.fwrk.VSEM;
quiverscale = 0.5;

%% The initial vector: 6 + 36 variables
y0 = zeros(42,1);
% STM concatenation after the 6-dim state
y0 = matrixToVector(y0, cst.orbit.STM0, 6, 6, 6);

%% Leg 1: EML2 orbit

% Initial conditions on a EML2 Lyapunov orbit
y0(1:6) = [-1.00171515716847 ; 0.00239911223400391   ; 0 ; -0.0308918369296315 ; -1.02460708873130  ;   0];
% Initial & final times in SE coordinates
t0 =  0.436887365666364;
tf = -0.0711211990619663;

%Integration
[tv, yv] = ode78_qbcp([t0 tf], y0, dcs, N1-1, inputType, outputType);

%Inverse the elements
teml2 = flipud(tv);
yeml2 = flipud(yv);


%% Leg 2: manifold leg

y0(1:6) = [-1.00171515717719 ; 0.002399112695912011  ; 0 ; -0.0308918461829167 ; -1.02460707457657  ;   0];
% Initial & final times in SE coordinates
t0 = 0.436887365666364;
tf = 2.96422997518981;

%Integration
[tman, yman] = ode78_qbcp([t0 tf], y0, dcs, N2-1, inputType, outputType);

%% Leg 3: SEML2 orbit

% Initial conditions on a SEML2 Lyapunov orbit
y0(1:6) = [-1.01028778775645 ; 0.00500510073123760   ; 0 ; -0.000626773575768337 ; -1.00384034451069  ;   0];
% Initial & final times in SE coordinates
t0 = 2.96422997518981;
tf = 4.48825566937480+1;

%Integration
[tseml2, yseml2] = ode78_qbcp([t0 tf], y0, dcs, N3-1, inputType, outputType);


%% Plot
figure(1)
hold on
grid on;
plot(yeml2(:,1),  yeml2(:,2),  'Color', 'k', 'LineWidth', 2);
plot(yman(:,1),   yman(:,2),   'Color', 'r', 'LineWidth', 2);
plot(yseml2(:,1), yseml2(:,2), 'Color', 'b', 'LineWidth', 2);
legend('EML2 orbit', 'Transfer leg', 'SEML2 orbit');

%% Initial and final values in NCSEM format
yeml2_0  = yman(1,:)';
yseml2_0 = yseml2(1,:)';

%% We change the scope: everything is in one type now
inputType  = cst.coord2.VNCSEM;
outputType = cst.coord2.VNCSEM;
dcs        = cst.fwrk.VSEM;

%% Back to the problem of the manifold
t0 = 0.436887365666364;
tf = 2.96422997518981;

%% Copy initial values
ymd = yman;
tmd = tman;

%--------------------------------------------------------------------------
% WARNING: for the last component, we integrate yseml2_0, since it is the
% target!
%--------------------------------------------------------------------------
ymd(end, 1:6) = yseml2(1,1:6);

%--------------------------------------------------------------------------
% WARNING: we need to update the STM at each patch point!
% Already done in ode78_qbcp, but you never know!
%--------------------------------------------------------------------------
for k = 1: size(ymd, 1)
    ymd(k,:) = matrixToVector(ymd(k,:), cst.orbit.STM0, 6, 6, 6);
end

%% Differential Correction

%--------------------------------------------------------------------------
% - tmdn, ymdn contains the new departure state & time
% - yma contains the arrival state
%--------------------------------------------------------------------------
isTimeFixed = true;
[tmdn, ymdn, yma] = diff_corr_level_I(tmd, ymd, N2, dcs, inputType, isTimeFixed);
%[tmdn, ymdn, yma] = diff_corr_level_II_with_constraints(tmd, ymd, N2, dcs, inputType, false);
%[tmdn, ymdn, yma] = diff_corr_level_II_fixed_tspan(tmd, ymd, N2, dcs, inputType);
%[tmdn, ymdn, yma] = diff_corr_level_II_fixed_tspan_mixed(tmd, ymd, N2, dcs, inputType);


%% Compute the cost (DV, DT)

cr3bp = init_CR3BP('SUN', 'EARTH_AND_MOON', default);
gamma_i = 0.010078240625297; %force a given gamma value, from C++ code

%--------------------------------------------------------------------------
% Compute DV in m/s
%--------------------------------------------------------------------------
DV = zeros(N2, 1);
for k = 1:N2
    DV(k) = 1e3*cr3bp.L*2*pi/cr3bp.T*gamma_i*norm(yma(k, 4:6) - ymdn(k, 4:6));
end

%--------------------------------------------------------------------------
% Compute DT in hours
%--------------------------------------------------------------------------
DT = zeros(N2, 1);
for k = 1:N2
    DT(k) = cr3bp.T/(2*pi*3600)*norm(tmdn(k) - tmd(k));
end