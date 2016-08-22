init;

%--------------------------------------------------------------------------
% Test of FMINCON for QBCP
%--------------------------------------------------------------------------
% Number of points on each legs
N1 = 100;
N2 = 2;
N3 = 100;

inputType  = cst.coord2.VNCSEM;
outputType = cst.coord2.VNCSEM;
dcs        = cst.fwrk.VSEM;
quiverscale = 0.5;

%% The initial vector: 6 + 36 variables
y0 = zeros(42,1);
% STM concatenation after the 6-dim state
y0 = matrixToVector(y0, cst.orbit.STM0, 6, 6, 6);


%% Leg 2: manifold leg
y0(1:6) = [ -7.909020665016300e-01;
            -2.057153010805557e-01;
             0.000000000000000e+00;
             2.397145759825049e+00;
             2.699054245560647e+00;
             0.000000000000000e+00];
% Initial & final times in SE coordinates
t0 = 4.500955883493010e-01;
tf = 2.901236913163497e+00;


%Integration
[tman, yman] = ode78_qbcp([t0 tf], y0, dcs, N2-1, inputType, outputType);

%% Target
yproj =  [4.311933504764909e-02;
         -4.935876996264280e-01;
          0.000000000000000e+00;
         -4.525820233983336e-01;
         -7.947431930418852e-01;
          0.000000000000000e+00];

%% Plot
figure(1)
hold on
grid on;
plot(yman(:,1),   yman(:,2),   'Color', 'r', 'LineWidth', 2);
legend('EML2 orbit', 'Transfer leg', 'SEML2 orbit');


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
ymd(end, 1:6) = yproj(1:6);

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
isTimeFixed = false;
%[tmdn, ymdn, yma] = diff_corr_level_II(tmd, ymd, N2, dcs, inputType, isTimeFixed);
[tmdn, ymdn, yma] = diff_corr_level_I(tmd, ymd, N2, dcs, inputType, isTimeFixed);


%[tmdn, ymdn, yma] = diff_corr_level_II_with_constraints(tmd, ymd, N2, dcs, inputType, false);
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
    ymdn(k, 4:6)
    yma(k, 4:6)
end

%--------------------------------------------------------------------------
% Compute DT in hours
%--------------------------------------------------------------------------
DT = zeros(N2, 1);
for k = 1:N2
    DT(k) = cr3bp.T/(2*pi*3600)*norm(tmdn(k) - tmd(k));
end