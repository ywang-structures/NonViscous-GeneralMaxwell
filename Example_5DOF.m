%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script demonstrates the procedure for obtaining the transformation
% matrix between the pair of 5-DOF systems. The configurations are devised
% to be more generic, including multiple connections with the ground.
% Dynamic response is simulated using the 1940 El Centro NS ground motion.
% The results of Section 8 of the manuscript are generated using this
% script.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear;

%% Define basic parameters
% Number of DOFs
N = 5;

% Mass
m = 10; % kg
M = diag([m, m, m, m, m]);

% Stiffness
k = 500; % N/m
K = [ 2*k  -k    0    0    0;
     -k    2*k  -k    0    0;
      0   -k    2*k  -k    0;
      0    0   -k    2*k  -k;
      0    0    0   -k    2*k ];

% Damping
n_damp = 2;
mu = [5, 10]; % 1/s
ct_I = 1; % Ns/m
ct_II = 1; % Ns/m

%% Generalized Maxwell model
% Maxwell damping matrix
G_gm{1} = [ ct_I   0     0;
             0    ct_I  -ct_I;
             0     0    ct_I;
             0     0     0;
             0     0     0 ];

G_gm{2} = [  0    0;
             0    0;
            -ct_II  0;
             ct_II  0;
             0   ct_II ];

% Decomposition G_gm = R_gm * D_gm
D_gm{1} = diag([ct_I, ct_I, ct_I]);
D_gm{2} = diag([ct_II, ct_II]);

R_gm{1} = [ 1  0   0;
            0  1  -1;
            0  0   1;
            0  0   0;
            0  0   0 ];

R_gm{2} = [  0  0;
             0  0;
            -1  0;
             1  0;
             0  1 ];

% Non-viscous damping matrix
Ct{1} = G_gm{1} * R_gm{1}';
Ct{2} = G_gm{2} * R_gm{2}';
rank_Ct = [rank(Ct{1}), rank(Ct{2})];

% Construct system matrices
% A_gm * xd_gm + B_gm * x_gm = F
n_x = 2*N + sum(rank_Ct);
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
  d_idx_begin(i) = 2*N + cum_rank_Ct(i) + 1;
  d_idx_end(i) = 2*N + cum_rank_Ct(i+1);
  
  A_gm(1:N, 1:N) = A_gm(1:N, 1:N) + Ct{i};
  A_gm(1:N, d_idx_begin(i):d_idx_end(i)) = -G_gm{i} / mu(i);
  A_gm(d_idx_begin(i):d_idx_end(i), 1:N) = -G_gm{i}' / mu(i);
  A_gm(d_idx_begin(i):d_idx_end(i), d_idx_begin(i):d_idx_end(i)) = D_gm{i} / mu(i)^2;
  B_gm(d_idx_begin(i):d_idx_end(i), d_idx_begin(i):d_idx_end(i)) = D_gm{i} / mu(i);
end

% Define input
groundAcc = load('ElCentro_NS.txt'); % g
g2m = 9.8;
fs = 200;
dt = 1/fs;
groundAcc = groundAcc' * g2m;
groundAcc = groundAcc(1:40*fs);
time = 0:dt:length(groundAcc)*dt-dt;
n_time = length(time);
p = -M * ones(N, 1) * groundAcc;
F = zeros(n_x, n_time);
F(1:N,:) = p; % assign input

% Solve state-space equation
% Convert from
%   A_gm * xd_gm + B_gm * x_gm = F
% To
%   xd_gm = Ass_gm * x_gm + Bss_gm * F
Ass_gm = -A_gm \ B_gm;
Bss_gm = inv(A_gm);
sys_gm = ss(Ass_gm, Bss_gm, eye(size(Ass_gm)), 0);
[x_gm,time] = lsim(sys_gm, F, time);
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
A_nv(1:N,N+1:2*N) = M;
A_nv(N+1:2*N,1:N) = M;
B_nv = zeros(n_x, n_x);
B_nv(1:N,1:N) = K;
B_nv(N+1:2*N,N+1:2*N) = -M;

for i = 1:n_damp
  % Eigenvalue decomposition of non-viscous damping matrix
  [R,D] = eigs(Ct{i}, rank_Ct(i));
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
[x_nv,time] = lsim(sys_nv, F, time);
x_nv = x_nv';

% Displacement
q_nv = x_nv(1:N, :);
% Velocity
qd_nv = x_nv(N+1:2*N, :);
% Condensed internal velocity
ud_nv_I = x_nv(d_idx_begin(1):d_idx_end(1), :);
ud_nv_II = x_nv(d_idx_begin(2):d_idx_end(2), :);

%% Transform matrix
for i = 1:n_damp
    T{i} = G_nv{i} \ G_gm{i};
end

% Transform damper velocity of GM to condensed internal velocity of NV
ud_nv_from_gm_I = T{1} * ud_gm_I;
ud_nv_from_gm_II = T{2} * ud_gm_II;


%% Check responses between NV and GM
fig_w = 800;
fig_h = 200;
figure('position', [100, 100, fig_w, fig_h])
dof = 1;
plot(time, q_nv(dof,:),'r'); hold on
plot(time, q_gm(dof,:),'b--');
legend('Displacement of NV','Displacement of GM')
xlabel('Time (sec)')
ylabel('Displacement (m)')

figure('position', [100, 100, fig_w, fig_h])
dof = 1;
plot(time, ud_nv_II(dof,:),'r'); hold on
plot(time, ud_nv_from_gm_II(dof,:),'b--');
legend('Condensed internal velocity of NV','Condensed internal velocity of NV obtained from damper velocity of GM')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')


%% Figure 10
% Damper velocity of GM
% Figure 10(a)
lw = 1.4;
lb = -4;
ub = 4;
figure('position', [100, 100, fig_w, fig_h])
plot(time, ud_gm_I(1,:),'r-','LineWidth',lw); hold on
plot(time, ud_gm_I(2,:),'b-.','LineWidth',lw); hold on
plot(time, ud_gm_I(3,:),'k:','LineWidth',lw); hold on
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
ylim([lb ub])

% Figure 10(b)
figure('position', [100, 100, fig_w, fig_h])
plot(time, ud_gm_II(1,:),'r','LineWidth',lw); hold on
plot(time, ud_gm_II(2,:),'b-.','LineWidth',lw);
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
ylim([lb ub])


%% Figure 11
% Condensed internal velocity of NV obtained by transformed damper velocity of GM
% Figure 11(a)
figure('position', [100, 100, fig_w, fig_h])
plot(time, ud_nv_from_gm_I(1,:),'r-','LineWidth',lw); hold on
plot(time, ud_nv_from_gm_I(2,:),'b-.','LineWidth',lw); hold on
plot(time, ud_nv_from_gm_I(3,:),'k:','LineWidth',lw); hold on
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
ylim([lb ub])

% Figure 11(b)
figure('position', [100, 100, fig_w, fig_h])
plot(time, ud_nv_from_gm_II(1,:),'r','LineWidth',lw); hold on
plot(time, ud_nv_from_gm_II(2,:),'b-.','LineWidth',lw);
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
ylim([lb ub])