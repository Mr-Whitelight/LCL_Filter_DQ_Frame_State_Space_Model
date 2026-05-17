% =========================================================================
% COMPLETE SCRIPT: Open-loop LCL + Full state feedback closed-loop
% =========================================================================
clear; clc; close all;

%% 1. SYMBOLIC DERIVATION OF OPEN-LOOP MODEL
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
A_sym_cl = jacobian(f, x);
B_sym_cl = jacobian(f, u);
y = [igd; igq];
C_sym_cl = jacobian(y, x);
D_sym_cl = jacobian(y, u);

%% 2. NUMERIC PARAMETERS (Table I)
L1v = 0.8e-3;  R1v = 50e-3;
L2v = 0.4e-3;  R2v = 50e-3;
Cv  = 30e-6;
wv  = 2*pi*50;

% Numeric base controller gains
kP_val  = 0.12;
kI_val  = 20.0;
k_i_val = 50.00;
k_v_val = 5.0;
k_ig_val = 10.0;   
fprintf('Controller gains: kP=%.3f, kI=%.1f, k_i=%.2f, k_v=%.2f, k_ig=%.2f\n', kP_val, kI_val, k_i_val, k_v_val, k_ig_val);


A_num = double(subs(A_sym_cl, [L1,R1,L2,R2,C,omega], [L1v,R1v,L2v,R2v,Cv,wv]));
B_num = double(subs(B_sym_cl, [L1,R1,L2,R2,C,omega], [L1v,R1v,L2v,R2v,Cv,wv]));
C_num = double(C_sym_cl);
D_num = double(D_sym_cl);

sys_ol = ss(A_num, B_num, C_num, D_num);

%% 3. OPEN-LOOP EIGENVALUE PLOT (with LCL parameters in title)
eig_open = eig(A_num);
figure(1);
plot(real(eig_open), imag(eig_open), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on; axis equal;
xlabel('Real part'); ylabel('Imaginary part');
title(sprintf(['Open-loop LCL filter eigenvalues ' ...
    '(L1=%.1f mH, L2=%.1f mH, C=%.1f µF, R1=%.3f Ω, R2=%.3f Ω)'], ...
    L1v*1e3, L2v*1e3, Cv*1e6, R1v, R2v));
line([0 0], ylim, 'Color','k','LineStyle','--');

fprintf('\n========== Open-loop eigenvalue analysis ==========\n');
print_eigenvalue_modes(eig_open, 'Open-loop');


%% SYMBOLIC DERIVATION OF CLOSED-LOOP MODEL (FULL STATE FEEDBACK)
% Parameters
syms L1_cl R1_cl L2_cl R2_cl C_cl omega_cl real
% Gains
syms kP kI k_i k_v k_ig real
% Plant states (dq)
syms id_cl iq_cl vCd_cl vCq_cl igd_cl igq_cl real
% Integrator states
syms Xd_cl Xq_cl real
% Reference inputs
syms igd_star_cl igq_star_cl real
% Disturbance inputs (grid voltages)
syms vgd_cl vgq_cl real

x_cl = [id_cl; iq_cl; vCd_cl; vCq_cl; igd_cl; igq_cl; Xd_cl; Xq_cl];
u_cl = [igd_star_cl; igq_star_cl];

% ----- Integrator dynamics -----
dXd_dt = igd_star_cl - igd_cl;
dXq_dt = igq_star_cl - igq_cl;

vd_cl = kP*(igd_star_cl - igd_cl) + kI*Xd_cl ...
        - k_i*id_cl - k_v*vCd_cl - k_ig*igd_cl;
vq_cl = kP*(igq_star_cl - igq_cl) + kI*Xq_cl ...
        - k_i*iq_cl - k_v*vCq_cl - k_ig*igq_cl;

% ----- Plant derivative expressions (with substituted vd_cl, vq_cl) -----
did_dt_cl   = (vd_cl - R1_cl*id_cl - vCd_cl + omega_cl*L1_cl*iq_cl) / L1_cl;
diq_dt_cl   = (vq_cl - R1_cl*iq_cl - vCq_cl - omega_cl*L1_cl*id_cl) / L1_cl;
dvCd_dt_cl  = (id_cl - igd_cl + omega_cl*C_cl*vCq_cl) / C_cl;
dvCq_dt_cl  = (iq_cl - igq_cl - omega_cl*C_cl*vCd_cl) / C_cl;
digd_dt_cl  = (vCd_cl - R2_cl*igd_cl - vgd_cl + omega_cl*L2_cl*igq_cl) / L2_cl;
digq_dt_cl  = (vCq_cl - R2_cl*igq_cl - vgq_cl - omega_cl*L2_cl*igd_cl) / L2_cl;


% ----- Full closed-loop derivative vector (8 states) -----
f_cl = [did_dt_cl; diq_dt_cl; dvCd_dt_cl; dvCq_dt_cl; ...
        digd_dt_cl; digq_dt_cl; dXd_dt; dXq_dt];

% Jacobians
Acl_sym = jacobian(f_cl, x_cl);
Bcl_sym = jacobian(f_cl, u_cl);
y_cl = [igd_cl; igq_cl];
Ccl_sym = jacobian(y_cl, x_cl);
Dcl_sym = jacobian(y_cl, u_cl);

Acl = double(subs(Acl_sym, ...
    [L1_cl, R1_cl, L2_cl, R2_cl, C_cl, omega_cl,  kP,  kI,  k_i,  k_v,  k_ig], ...
    [L1v,   R1v,   L2v,   R2v,   Cv,   wv,  kP_val, kI_val, k_i_val, k_v_val, k_ig_val]));
Bcl = double(subs(Bcl_sym, ...
    [L1_cl, R1_cl, L2_cl, R2_cl, C_cl, omega_cl,  kP,  kI,  k_i,  k_v,  k_ig], ...
    [L1v,   R1v,   L2v,   R2v,   Cv,   wv,  kP_val, kI_val, k_i_val, k_v_val, k_ig_val]));

Ccl = double(Ccl_sym);
Dcl = double(Dcl_sym);

sys_cl = ss(Acl, Bcl, Ccl, Dcl);

%% 5. CLOSED-LOOP EIGENVALUE PLOT (with controller gains in title)
eig_cl = eig(Acl);
figure(2);
plot(real(eig_cl), imag(eig_cl), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
grid on; axis equal;
xlabel('Real part'); ylabel('Imaginary part');
title(sprintf(['Closed-loop eigenvalues ' ...
    '(k_P=%.3f, k_I=%.1f, k_i=%.3f, k_v=%.3f, k_{ig}=%.3f)'], ...
    kP_val, kI_val, k_i_val, k_v_val, k_ig_val));
line([0 0], ylim, 'Color','k','LineStyle','--');

fprintf('\n========== Closed-loop eigenvalue analysis ==========\n');
print_eigenvalue_modes(eig_cl, 'Closed-loop');

%% 6. STEP RESPONSE PLOTS
figure(3);
step(sys_cl(1,1));
title(sprintf('Step response: i_{gd}^* -> i_{gd} (gains: k_P=%.3f, k_I=%.1f)', kP_val, kI_val));
grid on;



figure(4);
step(sys_cl(2,2));
title(sprintf('Step response: i_{gq}^* -> i_{gq} (gains: k_P=%.3f, k_I=%.1f)', kP_val, kI_val));
grid on;

%% 8. EIGENVALUE TRAJECTORIES vs. PROPORTIONAL GAIN k_P (kI fixed at base value)
kI_fixed = kI_val;   % use the base value from Section 4
kP_values = linspace(1, 600, 200);
eig_traj_kP = zeros(8, length(kP_values));

for idx = 1:length(kP_values)
    A_temp = Acl;  % Acl built with base kP_val, kI, k_i, k_v, k_ig
    kP_temp = kP_values(idx);
    % Update only the kP_val-dependent entries, keep all other gains as base
    A_temp(1,5) = -(kP_temp + k_ig_val) / L1v;   % igd term in d-axis
    A_temp(2,6) = -(kP_temp + k_ig_val) / L1v;   % igq term in q-axis
    % kI entries remain as base (already set in Acl)
    eig_traj_kP(:, idx) = eig(A_temp);
end





figure(6);
hold on;
for k = 1:8
    plot(real(eig_traj_kP(k,:)), imag(eig_traj_kP(k,:)), '.-');
end
grid on; axis equal;
xlabel('Real part'); ylabel('Imaginary part');
title(sprintf('Eigenvalue trajectories vs k_P (k_I = %.1f)', kI_fixed));
line([0 0], ylim, 'Color','k','LineStyle','--');
legend('Location','best');

% Find critical kP_val
max_real_kP = max(real(eig_traj_kP), [], 1);
crit_idx_kP = find(max_real_kP >= 0, 1, 'first');
if ~isempty(crit_idx_kP)
    fprintf('Critical k_P (first RHP crossing) = %.4f\n', kP_values(crit_idx_kP));
else
    disp('No RHP crossing for k_P sweep.');
end

%% 9. EIGENVALUE TRAJECTORIES vs. INTEGRAL GAIN k_I (kP_val fixed at base value)
kP_fixed = kP_val;   % use the base value
kI_values = linspace(1, 300000, 200);
eig_traj_kI = zeros(8, length(kI_values));

for idx = 1:length(kI_values)
    A_temp = Acl;
    kI_temp = kI_values(idx);
    % Update kI-dependent entries (integrator paths)
    A_temp(1,7) = kI_temp / L1v;   % Xd coefficient in d-axis
    A_temp(2,8) = kI_temp / L1v;   % Xq coefficient in q-axis
    % Ensure kP_val entries stay at base value (already set, but just to be safe)
    A_temp(1,5) = -(kP_fixed + k_ig_val) / L1v;
    A_temp(2,6) = -(kP_fixed + k_ig_val) / L1v;


    eig_traj_kI(:, idx) = eig(A_temp);
end

figure(7);
hold on;
for k = 1:8
    plot(real(eig_traj_kI(k,:)), imag(eig_traj_kI(k,:)), '.-');
end
grid on; axis equal;
xlabel('Real part'); ylabel('Imaginary part');
title(sprintf('Eigenvalue trajectories vs k_I (k_P = %.2f)', kP_fixed));
line([0 0], ylim, 'Color','k','LineStyle','--');
legend('Location','best');

% Find critical kI
max_real_kI = max(real(eig_traj_kI), [], 1);
crit_idx_kI = find(max_real_kI >= 0, 1, 'first');
if ~isempty(crit_idx_kI)
    fprintf('Critical k_I (first RHP crossing) = %.1f\n', kI_values(crit_idx_kI));
else
    disp('No RHP crossing for k_I sweep.');
end

disp('=== Analysis complete ===')

% =========================================================================
% Helper function to print eigenvalues with frequency, damping, and labels
% =========================================================================
function print_eigenvalue_modes(eigvals, sysName)
    format short g
    n = length(eigvals);
    % Sort by imaginary part (absolute) for grouping
    [~, idx] = sort(abs(imag(eigvals)), 'descend');
    eigvals = eigvals(idx);
    
    paired = false(n,1);
    modeNum = 0;
    fprintf('\n%s eigenvalues (rad/s) and their characteristics:\n', sysName);
    fprintf('%-8s %-20s %-12s %-12s %s\n', ...
        'Mode', 'Poles (rad/s)', 'Freq (Hz)', 'Damping', 'Interpretation');
    fprintf('%s\n', repmat('-', 1, 90));
    
    for k = 1:n
        if paired(k), continue; end
        p = eigvals(k);
        if imag(p) > 1e-6   % complex pole with positive imag
            % find its conjugate
            conj_idx = find(abs(eigvals - conj(p)) < 1e-6);
            if ~isempty(conj_idx)
                paired([k, conj_idx]) = true;
                w = abs(imag(p));
                freqHz = w/(2*pi);
                sigma = real(p);
                zeta = -sigma/abs(p);
                % Label based on frequency
                if freqHz > 1000
                    label = 'LCL resonance (high freq)';
                elseif freqHz > 40 && freqHz < 60
                    label = 'Grid-frequency mode (50 Hz)';
                elseif freqHz < 1
                    label = 'Integrator mode (slow)';
                else
                    label = 'Other';
                end
                modeNum = modeNum + 1;
                fprintf('%2d   %8.2f ± %8.2fi  %10.2f  %10.4f  %s\n', ...
                    modeNum, sigma, w, freqHz, zeta, label);
            end
        else   % real pole
            paired(k) = true;
            freqHz = 0;
            sigma = real(p);
            zeta = 1; % for real pole, damping = 1
            if abs(sigma) < 10
                label = 'Integrator (real)';
            else
                label = 'Real pole';
            end
            modeNum = modeNum + 1;
            fprintf('%2d   %8.2f               %10.2f  %10.4f  %s\n', ...
                modeNum, sigma, freqHz, zeta, label);
        end
    end
end
