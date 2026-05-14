% =========================================================================
% COMPLETE SCRIPT: Open-loop LCL + Full state feedback closed-loop
% =========================================================================
clear; clc;

%% 1. SYMBOLIC DERIVATION OF OPEN-LOOP MODEL (transformerless)
syms L1 R1 L2 R2 C omega real
syms id iq vCd vCq igd igq real
syms vd vq vgd vgq real

x = [id; iq; vCd; vCq; igd; igq];
u = [vd; vq; vgd; vgq];

% Derivative expressions
did_dt   = (vd - R1*id - vCd + omega*L1*iq) / L1;
diq_dt   = (vq - R1*iq - vCq - omega*L1*id) / L1;
dvCd_dt  = (id - igd + omega*C*vCq) / C;
dvCq_dt  = (iq - igq - omega*C*vCd) / C;
digd_dt  = (vCd - R2*igd - vgd + omega*L2*igq) / L2;
digq_dt  = (vCq - R2*igq - vgq - omega*L2*igd) / L2;

f = [did_dt; diq_dt; dvCd_dt; dvCq_dt; digd_dt; digq_dt];

% Jacobians
A_sym = jacobian(f, x);
B_sym = jacobian(f, u);
y = [igd; igq];
C_sym = jacobian(y, x);
D_sym = jacobian(y, u);

%% 2. NUMERIC PARAMETERS (Table I)
L1v = 1e-3;  R1v = 50e-3;
L2v = 400e-6;   R2v = 50e-3;
Cv  = 32.6e-6;
wv  = 2*pi*50;

A_num = double(subs(A_sym, [L1,R1,L2,R2,C,omega], [L1v,R1v,L2v,R2v,Cv,wv]));
B_num = double(subs(B_sym, [L1,R1,L2,R2,C,omega], [L1v,R1v,L2v,R2v,Cv,wv]));
C_num = double(C_sym);
D_num = double(D_sym);

%% 3. OPEN-LOOP EIGENVALUE PLOT
eig_open = eig(A_num);
figure(1);
plot(real(eig_open), imag(eig_open), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on; axis equal;
xlabel('Real part'); ylabel('Imaginary part');
title('Open-loop LCL filter eigenvalues (6 states)');
line([0 0], ylim, 'Color','k','LineStyle','--');
disp('Open-loop eigenvalues:');
disp(eig_open);

%% 4. FULL STATE FEEDBACK CONTROLLER DESIGN
% PI + state feedback gains
kP  = 0.5;
kI  = 1.0;
k_i = 0;   % converter-side current feedback
k_v = 0;   % capacitor voltage feedback
k_ig = 0;  % grid-side current feedback

% Build augmented closed-loop A matrix (8x8)
Acl = zeros(8,8);
Acl(1:6,1:6) = A_num;   % original plant

% d-axis (row 1)
Acl(1,1) = -(R1v + k_i) / L1v;
Acl(1,2) =  wv;
Acl(1,3) = -(1 + k_v) / L1v;
Acl(1,5) = -(kP + k_ig) / L1v;
Acl(1,7) =  kI / L1v;

% q-axis (row 2)
Acl(2,1) = -wv;
Acl(2,2) = -(R1v + k_i) / L1v;
Acl(2,4) = -(1 + k_v) / L1v;
Acl(2,6) = -(kP + k_ig) / L1v;
Acl(2,8) =  kI / L1v;

% Rows 3-6: unchanged
% (already set via Acl(1:6,1:6) = A_num;)

% Integrator rows
Acl(7,5) = -1;   % Xd derivative
Acl(8,6) = -1;   % Xq derivative

% B matrix for references (igd*, igq*)
Bcl = zeros(8,2);
Bcl(1,1) = kP/L1v;
Bcl(2,2) = kP/L1v;
Bcl(7,1) = 1;
Bcl(8,2) = 1;

% Output matrix
Ccl = [0 0 0 0 1 0 0 0;
       0 0 0 0 0 1 0 0];
Dcl = zeros(2,2);

sys_cl = ss(Acl, Bcl, Ccl, Dcl);

%% 5. CLOSED-LOOP EIGENVALUE PLOT
eig_cl = eig(Acl);
figure(2);
plot(real(eig_cl), imag(eig_cl), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
grid on; axis equal;
xlabel('Real part'); ylabel('Imaginary part');
title('Closed-loop eigenvalues (full state feedback, 8 states)');
line([0 0], ylim, 'Color','k','LineStyle','--');
disp('Closed-loop eigenvalues:');
disp(eig_cl);

%% 6. STEP RESPONSE PLOTS
figure(3);
step(sys_cl(1,1));
title('Step response: i_{gd}^* -> i_{gd}');
grid on;

figure(4);
step(sys_cl(2,2));
title('Step response: i_{gq}^* -> i_{gq}');
grid on;

%% 7. EIGENVALUE TRAJECTORIES vs. CAPACITOR VOLTAGE FEEDBACK k_v
k_v_values = linspace(0, 0.2, 30);
eig_traj = zeros(8, length(k_v_values));

for idx = 1:length(k_v_values)
    A_temp = Acl;
    k_v_temp = k_v_values(idx);
    A_temp(1,3) = -(1 + k_v_temp) / L1v;
    A_temp(2,4) = -(1 + k_v_temp) / L1v;
    eig_traj(:, idx) = eig(A_temp);
end

figure(5);
hold on;
for k = 1:8
    plot(real(eig_traj(k,:)), imag(eig_traj(k,:)), '.-');
end
grid on; axis equal;
xlabel('Real part'); ylabel('Imaginary part');
title('Eigenvalue trajectories as k_v varies (active damping)');
line([0 0], ylim, 'Color','k','LineStyle','--');
legend('Location','best');

disp('=== Analysis complete ===')

