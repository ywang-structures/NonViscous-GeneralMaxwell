%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using a pair of 3-DOF systems, this script first demonstrates the
% procedure for obtaining the transformation matrix that relates the
% internal velocities of exponential non-viscous damping to the damper
% velocities of the generalized Maxwell model. Additionally, it analyzes
% the modal response of both oscillatory and non-oscillatory modes in
% exponential non-viscous damping. The results presented in Section 7 of
% the manuscript are generated using this script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear;

%% Define basic parameters
% Number of DOFs
N = 3;

% Mass
m1 = 5; % kg
m2 = 5;
m3 = 5;
M = diag([m1, m2, m3]);

% Stiffness
k1 = 800; % N/m
k2 = 1000;
k3 = 1200;
K = [k1 + k2, -k2, 0;
     -k2, k2 + k3, -k3;
      0, -k3, k3];

% Damping
n_damp = 2;
mu = [5, 10]; % 1/s
ct_I = [1, 2]; % Ns/m
ct_II = [5, 10];

%% Generalized Maxwell model
% Maxwell damping matrix
G_gm{1} = [ct_I(1), -ct_I(2);
            0, ct_I(2);
            0, 0];
G_gm{2} = [-ct_II(1), 0;
           ct_II(1), -ct_II(2);
           0, ct_II(2)];

% Decomposition G_gm = R_gm * D_gm
D_gm{1} = diag([ct_I(1), ct_I(2)]);
D_gm{2} = diag([ct_II(1), ct_II(2)]);
R_gm{1} = [1, -1;
            0,  1;
            0,  0];
R_gm{2} = [-1, 0;
            1, -1;
            0,  1];

% Non-viscous damping matrix
Ct{1} = G_gm{1} * R_gm{1}';
Ct{2} = G_gm{2} * R_gm{2}';
rank_Ct = [rank(Ct{1}), rank(Ct{2})];

% Construct system matrices
% A_gm * xd_gm + B_gm * x_gm = F
n_x = 2 * N + sum(rank_Ct);
A_gm = zeros(n_x, n_x);
A_gm(1:N, N+1:2*N) = M;
A_gm(N+1:2*N, 1:N) = M;
B_gm = zeros(n_x, n_x);
B_gm(1:N, 1:N) = K;
B_gm(N+1:2*N, N+1:2*N) = -M;

% Define index for convenience
d_idx_begin = zeros(n_damp, 1);
d_idx_end = zeros(n_damp, 1);
cum_rank_Ct = [0, cumsum(rank_Ct)];

% Loop through dampers
for i = 1:n_damp
    d_idx_begin(i) = 2 * N + cum_rank_Ct(i) + 1;
    d_idx_end(i) = 2 * N + cum_rank_Ct(i+1);
    
    A_gm(1:N, 1:N) = A_gm(1:N, 1:N) + Ct{i};
    A_gm(1:N, d_idx_begin(i):d_idx_end(i)) = -G_gm{i} / mu(i);
    A_gm(d_idx_begin(i):d_idx_end(i), 1:N) = -G_gm{i}' / mu(i);
    A_gm(d_idx_begin(i):d_idx_end(i), d_idx_begin(i):d_idx_end(i)) = D_gm{i} / mu(i)^2;
    B_gm(d_idx_begin(i):d_idx_end(i), d_idx_begin(i):d_idx_end(i)) = D_gm{i} / mu(i);
end

% Define initial condition and time duration
x0 = [ones(N,1); zeros(n_x-N,1)];
duration = 20; % sec
fs = 200; % Hz
dt = 1/fs;
time = 0:dt:duration-dt;
n_time = length(time);
F = zeros(n_x, n_time);

% Solve state-space equation
% Convert from:
%   A_gm * xd_gm + B_gm * x_gm = F
% To:
%   xd_gm = Ass_gm * x_gm + Bss_gm * F
Ass_gm = -A_gm \ B_gm;
Bss_gm = inv(A_gm);
sys_gm = ss(Ass_gm, Bss_gm, eye(size(Ass_gm)), 0);
[x_gm, time] = lsim(sys_gm, F, time, x0);
x_gm = x_gm';

% Displacement
q_gm = x_gm(1:N, :);
% Velocity
qd_gm = x_gm(N+1:2*N, :);
% Damper velocity
ud_gm_I = x_gm(d_idx_begin(1):d_idx_end(1), :);
ud_gm_II = x_gm(d_idx_begin(2):d_idx_end(2), :);

%% Exponential non-viscous damping
% Construct system matrices
A_nv = zeros(n_x, n_x);
A_nv(1:N, N+1:2*N) = M;
A_nv(N+1:2*N, 1:N) = M;
B_nv = zeros(n_x, n_x);
B_nv(1:N, 1:N) = K;
B_nv(N+1:2*N, N+1:2*N) = -M;

for i = 1:n_damp
    % Eigenvalue decomposition of non-viscous damping matrix
    [R, D] = eigs(Ct{i}, rank_Ct(i));
    R_nv{i} = R;
    D_nv{i} = D;
    G_nv{i} = R_nv{i} * D_nv{i};

    A_nv(1:N, 1:N) = A_nv(1:N, 1:N) + Ct{i};
    A_nv(1:N, d_idx_begin(i):d_idx_end(i)) = -G_nv{i} / mu(i);
    A_nv(d_idx_begin(i):d_idx_end(i), 1:N) = -G_nv{i}' / mu(i);
    A_nv(d_idx_begin(i):d_idx_end(i), d_idx_begin(i):d_idx_end(i)) = D_nv{i} / mu(i)^2;
    B_nv(d_idx_begin(i):d_idx_end(i), d_idx_begin(i):d_idx_end(i)) = D_nv{i} / mu(i);
end

% Solve state-space equation
Ass_nv = -A_nv \ B_nv;
Bss_nv = inv(A_nv);
sys_nv = ss(Ass_nv, Bss_nv, eye(size(Ass_nv)), 0);
[x_nv, time] = lsim(sys_nv, F, time, x0);
x_nv = x_nv';

% Displacement
q_nv = x_nv(1:N, :);
% Velocity
qd_nv = x_nv(N+1:2*N, :);
% Condensed internal velocity
ud_nv_I = x_nv(d_idx_begin(1):d_idx_end(1), :);
ud_nv_II = x_nv(d_idx_begin(2):d_idx_end(2), :);

%% Transformation matrix
for i = 1:n_damp
    T{i} = G_nv{i} \ G_gm{i};
end

% Transform damper velocity of GM to condensed internal velocity of NV
ud_nv_from_gm_I = T{1} * ud_gm_I;
ud_nv_from_gm_II = T{2} * ud_gm_II;

%% Figure 6 damper velocities and condensed internal velocities
fig_w = 800;
fig_h = 350;
lb = -5;
ub = 5;

% Figure 6(a)
figure('position', [100, 100, fig_w, fig_h])
dof = 1;
subplot(2,1,1)
plot(time, ud_gm_I(dof,:),'b'); hold on
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
ylim([lb ub])
subplot(2,1,2)
plot(time, ud_nv_I(dof,:),'k'); hold on
plot(time, ud_nv_from_gm_I(dof,:),'b:', 'LineWidth',1.5);
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
ylim([lb ub])

% Figure 6(b)
dof = 2;
figure('position', [100, 100, fig_w, fig_h])
subplot(2,1,1)
plot(time, ud_gm_II(dof,:),'b'); hold on
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
ylim([lb ub])
subplot(2,1,2)
plot(time, ud_nv_II(dof,:),'k'); hold on
plot(time, ud_nv_from_gm_II(dof,:),'b:', 'LineWidth', 2);
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
ylim([lb ub])

%% Eigenvalue and eigenvector analysis
% Table 2
[eigvec_gm, eigval_gm] = eigs(B_gm, -A_gm, n_x,'sm');
eigval_gm = diag(eigval_gm);
[~,idx_sort] = sort(abs(real(eigval_gm)),'ascend');
eigval_gm = eigval_gm(idx_sort);
eigvec_gm = eigvec_gm(:, idx_sort);
reorder = [1:2:2*N, 2:2:2*N, 2*N+1:n_x];
eigval_gm = eigval_gm(reorder);
eigvec_gm = eigvec_gm(:, reorder);

[eigvec_nv, eigval_nv] = eigs(B_nv, -A_nv, n_x, 'sm');
eigval_nv = diag(eigval_nv);
[~,idx_sort] = sort(abs(real(eigval_nv)),'ascend');
eigval_nv = eigval_nv(idx_sort);
eigvec_nv = eigvec_nv(:, idx_sort);
eigval_nv = eigval_nv(reorder);
eigvec_nv = eigvec_nv(:, reorder);

% Normalize by the largest entry of the first N entries
for i = 1:n_x
    [~, qth] = max(abs(eigvec_nv(1:N,i)));

    qth_val_nv = eigvec_nv(qth,i);
    eigvec_nv(:,i) = eigvec_nv(:,i) / qth_val_nv;

    qth_val_gm = eigvec_gm(qth,i);
    eigvec_gm(:,i) = eigvec_gm(:,i) / qth_val_gm;
end
H_from_eigvec = eigvec_nv / eigvec_gm; % Eq (108)

% Confirm that the transformation matrix obtained by the two approaches is
% the same
for i = 1:n_damp
    T_from_eigvec{i} = H_from_eigvec(d_idx_begin(i):d_idx_end(i), d_idx_begin(i):d_idx_end(i));
    is_T_same(i) = norm(T{i} - T_from_eigvec{i}) < 1e-9;
end

%% Section 7.2 Modal response of exponential non-viscous damping system
n_mode = n_x;
z_nv_mode = zeros(n_mode, n_time); % SDOF modal response
z0 = eigvec_nv \ x0; % initial condition for each mode
x_nv_ms = zeros(n_x, n_time); % response obtained from modal superposition
for r = 1:n_mode
    A_nv_mode = eigval_nv(r);
    B_nv_mode = transpose(eigvec_nv(:,r));
    C_nv_mode = eye(1);
    D_nv_mode = 0;
    sys = ss(A_nv_mode, B_nv_mode, C_nv_mode, D_nv_mode);
    z_nv_mode(r,:) = lsim(sys, F, time, z0(r));
    x_nv_mode = eigvec_nv(:,r) * z_nv_mode(r,:);

    % Modal response of displacement
    q_nv_mode(:,:,r) = x_nv_mode(1:N,:);
    % Modal response of condensed internal velocity
    ud_nv_I_mode(:,:,r) = x_nv_mode(d_idx_begin(1):d_idx_end(1),:);
    ud_nv_II_mode(:,:,r) = x_nv_mode(d_idx_begin(2):d_idx_end(2),:);

    % modal superposition
    x_nv_ms = x_nv_ms + x_nv_mode;
end
% Check the modal superposition response
is_x_same = norm(x_nv_ms - x_nv) < 1e-9;

% Displacement
q_nv_ms = x_nv_ms(1:N,:);
% Velocity
qd_nv_ms = x_nv_ms(N+1:2*N,:);
% Condensed internal velocity
ud_nv_ms_I = x_nv_ms(d_idx_begin(1):d_idx_end(1),:);
ud_nv_ms_II = x_nv_ms(d_idx_begin(2):d_idx_end(2),:);

%%  Figure 7 Modal response of displacement
fig_w = 800;
fig_h = 150;
warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
lw = 1;

dof = 3;
%  Figure 7(a) total response
figure('position', [100, 100, fig_w, fig_h])
plot(time, q_nv_ms(dof,:),'k','LineWidth',lw);hold on
xlabel('Time (sec)')
ylabel('Displacement (m)')
xlim([0 3])

%  Figure 7(b) oscillatory modes
figure('position', [100, 100, fig_w, fig_h])
plot(time, 2*real(q_nv_mode(dof,:,1)),'Color',[.5 .5 .5],'LineWidth',lw);hold on
plot(time, 2*real(q_nv_mode(dof,:,2)),'-r','LineWidth',lw);
plot(time, 2*real(q_nv_mode(dof,:,3)),'-b','LineWidth',lw);
xlabel('Time (sec)')
ylabel('Displacement (m)')
legend('r = 1','r = 2','r = 3')
xlim([0 3])

%  Figure 7(c) non-oscillatory modes
figure('position', [100, 100, fig_w, fig_h])
plot(time, q_nv_mode(dof,:,7),'-.','Color',[.5 .5 .5],'LineWidth',lw);hold on
plot(time, q_nv_mode(dof,:,8),'-.r','LineWidth',lw);
plot(time, q_nv_mode(dof,:,9),'-.b','LineWidth',lw);
plot(time, q_nv_mode(dof,:,10),'-.','Color',[0 .5 0],'LineWidth',lw);
xlabel('Time (sec)')
ylabel('Displacement (m)')
legend('r = 7','r = 8','r = 9','r = 10')
xlim([0 3])

%%  Figure 8 Modal response of condensed internal velocity
dof = 2;
%  Figure 8(a) total response
figure('position', [100, 100, fig_w, fig_h])
plot(time, ud_nv_ms_II(dof,:),'k','LineWidth',lw);hold on
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
xlim([0 3])

%  Figure 8(b) oscillatory modes
figure('position', [100, 100, fig_w, fig_h])
plot(time, 2*real(ud_nv_II_mode(dof,:,1)),'Color',[.5 .5 .5],'LineWidth',lw);hold on
plot(time, 2*real(ud_nv_II_mode(dof,:,2)),'-r','LineWidth',lw);
plot(time, 2*real(ud_nv_II_mode(dof,:,3)),'-b','LineWidth',lw);
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('r = 1','r = 2','r = 3')
xlim([0 3])

%  Figure 8(c) non-oscillatory modes
figure('position', [100, 100, fig_w, fig_h])
plot(time, ud_nv_II_mode(dof,:,7),'-.','Color',[.5 .5 .5],'LineWidth',lw);hold on
plot(time, ud_nv_II_mode(dof,:,8),'-.r','LineWidth',lw);
plot(time, ud_nv_II_mode(dof,:,9),'-.b','LineWidth',lw);
plot(time, ud_nv_II_mode(dof,:,10),'-.','Color',[0 .5 0],'LineWidth',lw);
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('r = 7','r = 8','r = 9','r = 10')
xlim([0 3])
ylim([-2 0.1])
