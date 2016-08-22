%--------------------------------------------------------------------------
% Test of FMINCON for QBCP
%--------------------------------------------------------------------------

% Number of points on each legs
N1 = 20;
N2 = 15;
N3 = 20;

%% Leg 1: EML2 orbit

% Initial conditions on a EML2 Lyapunov orbit
y0 = [-1.00171515716847 ; 0.00239911223400391   ; 0 ; -0.0308918369296315 ; -1.02460708873130  ;   0];
% Initial & final times in SE coordinates
t0 =  0.436887365666364;
tf = -0.0711211990619663;

%Integration
[tv, yv] = ode78_qbcp([t0 tf], y0, 0, N1-1);

%Inverse the elements
teml2 = flipud(tv);
yeml2 = flipud(yv);

%% Leg 2: manifold leg

y0 = [-1.00171515717719 ; 0.002399112695912011  ; 0 ; -0.0308918461829167 ; -1.02460707457657  ;   0];
% Initial & final times in SE coordinates
t0 = 0.436887365666364;
tf = 2.96422997518981;

%Integration
[tman, yman] = ode78_qbcp([t0 tf], y0, 0, N2-1);

%% Leg 3: SEML2 orbit

% Initial conditions on a SEML2 Lyapunov orbit
y0 = [-1.01028778775645 ; 0.00500510073123760   ; 0 ; -0.000626773575768337 ; -1.00384034451069  ;   0];
% Initial & final times in SE coordinates
t0 = 2.96422997518981;
tf = 4.48825566937480+1;

%Integration
[tseml2, yseml2] = ode78_qbcp([t0 tf], y0, 0, N3-1);


%% Plot
figure(1)
hold on
grid on;
plot(yeml2(:,1),  yeml2(:,2),  'Color', 'k', 'LineWidth', 2);
plot(yman(:,1),   yman(:,2),   'Color', 'r', 'LineWidth', 2);
plot(yseml2(:,1), yseml2(:,2), 'Color', 'b', 'LineWidth', 2);
legend('EML2 orbit', 'Transfer leg', 'SEML2 orbit');

%% Build the variable vector
% Build the state x0 = [y1,...,yN, t1,...tN]
x0 = [];
for i = 1:N2
    x0 = [x0 ; yman(i,:)'];
end
for i = 1:N2
    x0 = [x0 ; tman(i)];
end

%% Initial and final values
yeml2  = [-1.00171515717719 ; 0.002399112695912011  ; 0 ; -0.0308918461829167 ; -1.02460707457657  ;   0];
yseml2 = [-1.01028778775645 ; 0.00500510073123760   ; 0 ; -0.000626773575768337 ; -1.00384034451069  ;   0];
t0 = 0.436887365666364;
tf = 2.96422997518981;

%% Fmincon
options = optimoptions('fmincon','Display','iter','Algorithm', 'interior-point');
%options = optimoptions('fmincon','Display','iter','Algorithm', 'sqp');

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

nonlcon = @(x)nonlinconditions_qbcp(x, yeml2, yseml2, t0, tf, N2);
fun     = @(x)minconditions_qbcp(x, yeml2, yseml2, N2);
x       = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

%% Output
yvv2 = [];
for k = 1:N2
    % State at tk
    yvv2(k,:) = x((k-1)*6+1:k*6)';
end

figure(1)
hold on
plot(yvv2(:,1), yvv2(:,2), 'o');


%% Integration
y0 = yvv2(1,:);
[tman, yman] = ode78_qbcp([x(6*N2+1) x(7*N2)], y0, 0, 100);

figure(1)
hold on
plot(yman(:,1),   yman(:,2),   'Color', 'm', 'LineWidth', 2);