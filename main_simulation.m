% =========================================================================
% Stochastic ILC Simulation (Section 5)
% Comparison of Proposed Algorithm 3.1, SA-Shen [33], and SA-Cheng [10]

clear; clc; close all;

set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultAxesFontSize',10)
set(groot,'defaultAxesLabelFontSizeMultiplier',1) 
set(groot,'defaultFigurePosition',[600 500 400 300])
rng(42);

%% 1. System Parameters and Initialization
N = 51;         % Time duration
n = 4;          % State 
p = 3;          % Input 
q = 2;          % Output 
num_iter = 51;  % Iterations
Ks = 3000;        % Monte Carlo trials 

% Initialize 
A = cell(N, 1);
B = cell(N, 1);
C = cell(N, 1);

for t_idx = 1:N
    t = t_idx - 1; 
    A{t_idx} = [0.12, 0, 0, 0;
                0.03*exp(0.03*t), -0.2, -0.1, 0.03/(t+1);
                0, 0.05, 0, 0.03*cos(3*t);
                -0.03*t, 0, 0, -0.1];
            
    B{t_idx} = [0.3, 0, 0;
                0, 0.5, -0.3*t;
                sin(0.3*t), 0.1, 0.3;
                0, 5+2*cos(2*t), 5*t+3];
            
    t_C = t_idx; 
    C{t_idx} = [3, 0, 0.3*cos(t_C-1), 0;
                0.3*(t_C-1), 5, 0, 0.3];
end

% Prescribed Trajectory
yd = zeros(N*q, 1);
for t_idx = 1:N
    t = t_idx;
    yd((t_idx-1)*q+1 : t_idx*q) = [10*(t/50)^2 * (1 - t/50);
                                   1.5*sin(0.04*pi*t)];
end

%% 2. Lifted Matrices (H_sys = G, M = X, Gw = D in paper)
H_sys = zeros(N*q, N*p);
M = zeros(N*q, n);
Gw = zeros(N*q, N*n);

for i = 1:N
    val_M = C{i};
    for k = i-1:-1:1
        val_M = val_M * A{k};
    end
    M((i-1)*q+1:i*q, :) = val_M;
    
    for j = 1:i
        val_H = C{i};
        for k = i-1:-1:j
            val_H = val_H * A{k};
        end
        val_H = val_H * B{j};
        H_sys((i-1)*q+1:i*q, (j-1)*p+1:j*p) = val_H;
        
        val_Gw = C{i};
        for k = i-1:-1:j
            val_Gw = val_Gw * A{k};
        end
        Gw((i-1)*q+1:i*q, (j-1)*n+1:j*n) = val_Gw;
    end
end

%% 3. Compute (u_inf) 
u_inf = zeros(N*p, 1);
Gamma_inf = (H_sys'*H_sys + 0.15*eye(N*p)) \ H_sys'; 
for k = 1:200
    Y_det = H_sys * u_inf;
    u_inf = u_inf + Gamma_inf * (yd - Y_det);
end

%% 4. Algorithm Setup
% Weighting Matrices
Q_weight = 1*eye(N*q);
R_weight = 0.01 * eye(N*p);
Q_sqrt = sqrtm(Q_weight);
Q_inv_sqrt = inv(Q_sqrt);

% Iteration-independent matrix K (From Table 2)
K_mat = Q_sqrt * H_sys * inv(H_sys' * Q_weight * H_sys + R_weight) * H_sys' * Q_sqrt;

% Baseline matrix for Comparison Algorithms
L_base = zeros(N*p, N*q);
for i = 1:N
    CB = C{i} * B{i};
    L_base((i-1)*p+1:i*p, (i-1)*q+1:i*q) = pinv(CB);
end

% 1. SA-Shen Parameters
L_shen = 2 * L_base; 

% 2. SA-Cheng Parameters
theta_cheng = 1.5; 
alpha_cheng = 0.2;
N0_cheng = N * p;

% Tracking variables
err_out_prop = zeros(num_iter, 1); err_in_prop = zeros(num_iter, 1);
err_out_shen = zeros(num_iter, 1); err_in_shen = zeros(num_iter, 1);
err_out_cheng = zeros(num_iter, 1); err_in_cheng = zeros(num_iter, 1);

% Variables to store output for plotting
Y_prop_k5 = zeros(N*q, 1);
Y_prop_k50 = zeros(N*q, 1);

%% 5. Monte Carlo Simulation
fprintf('Starting Monte Carlo Simulation (%d trials)...\n', Ks);
for trial = 1:Ks
    % 5.1 Noise Variances 
    % p_w(t) and p_v(t) 
    pw_sqrt = 0.002 + 0.008 * rand(N, 1); 
    pv_sqrt = 0.002 + 0.008 * rand(N, 1);
    
    Pw = diag(kron(pw_sqrt.^2, ones(n, 1)));
    Pv = diag(kron(pv_sqrt.^2, ones(q, 1)));
    Px = (0.002^2) * eye(n);
    
    % Compute Disturbance Covariance Matrix F (Table 2)
    F_mat = Gw * Pw * Gw' + M * Px * M' + Pv;
    
    % Initial inputs
    u_prop = zeros(N*p, 1);
    u_shen = zeros(N*p, 1);
    u_cheng = zeros(N*p, 1);
    
    % Proposed Algorithm 3.1 Initializations
    % H_0 = Q^{1/2} * Y_d * Y_d^T * Q^{1/2} (Since U_0 = 0)
%     H_k = 10000*Q_sqrt * (yd * yd') * Q_sqrt;
    H_k = Q_sqrt * (yd * yd') * Q_sqrt;
    
    % SA-Shen variables
    s_shen = 1;
    E_prev_shen = zeros(N*q, 1);
    
    % SA-Cheng variables
    eta_cheng = 1;
    E_prev_cheng = zeros(N*q, 1);
    p_cheng = (1/N0_cheng) * ones(N0_cheng, 1);
    
    for k = 1:num_iter
        % Generate noise values
        W = randn(N*n, 1) .* kron(pw_sqrt, ones(n, 1));
        V = randn(N*q, 1) .* kron(pv_sqrt, ones(q, 1));
        x0_prime = 0.01 * randn(n, 1);
        
        % Lifted disturbance dk
        dk = M * x0_prime + Gw * W + V;
        
        % =========================================================
        % 5.2 Algorithm 3.1
        Y_prop = H_sys * u_prop + dk;
        E_prop = yd - Y_prop;
        
        err_out_prop(k) = err_out_prop(k) + norm(E_prop)^2;
        err_in_prop(k)  = err_in_prop(k)  + norm(u_inf - u_prop)^2;
        
        % Compute Gamma_k (Eq. 3.15)
        Inv_term = inv(Q_inv_sqrt * H_k * Q_inv_sqrt + F_mat);
        Gamma_k = inv(H_sys'*Q_weight*H_sys + R_weight) * H_sys' * Q_sqrt * H_k * Q_inv_sqrt * Inv_term;
        
        % Update Control Input (Eq. 3.7)
        u_prop = u_prop + Gamma_k * E_prop;
        
        % Compute L_k (Table 2)
        L_k = H_k * Q_inv_sqrt * Inv_term * Q_inv_sqrt * H_k;
        
        % Update H_{k+1} (Eq. 3.16)
        H_k = H_k - K_mat * L_k - L_k * K_mat + K_mat * L_k * K_mat;
        
        % Store outputs 
        if k == 5,  Y_prop_k5 = Y_prop;  end
        
        % =========================================================
        % 5.3 SA-Shen Algorithm
        % =========================================================
        Y_shen = H_sys * u_shen + dk;
        E_shen = yd - Y_shen;
        
        err_out_shen(k) = err_out_shen(k) + norm(E_shen)^2;
        err_in_shen(k)  = err_in_shen(k)  + norm(u_inf - u_shen)^2;
        
        if k > 1 && (E_shen' * E_prev_shen < 0)
            s_shen = s_shen + 1;
        end
        u_shen = u_shen + (1/s_shen) * L_shen * E_shen;
        E_prev_shen = E_shen;
        
        % =========================================================
        % 5.4 SA-Cheng Algorithm
        % =========================================================
        Y_cheng = H_sys * u_cheng + dk;
        E_cheng = yd - Y_cheng;
        
        err_out_cheng(k) = err_out_cheng(k) + norm(E_cheng)^2;
        err_in_cheng(k)  = err_in_cheng(k)  + norm(u_inf - u_cheng)^2;
        
        if k > 1 && (E_cheng' * E_prev_cheng < 0)
            eta_cheng = eta_cheng + 1;
        end
        
        grad_k = L_base * E_cheng; 
        
        if k == 1
            p_cheng = (1/N0_cheng) * ones(N0_cheng, 1);
        else
            grad_prev = L_base * E_prev_cheng;
            sgn_match = (grad_prev .* grad_k) > 0;
            m_k = sum(sgn_match);
            
            lambda_k = ones(N0_cheng, 1) / N0_cheng;
            if m_k > 0
                lambda_k(sgn_match) = 1 / m_k;
            end
            p_cheng = alpha_cheng * p_cheng + (1 - alpha_cheng) * lambda_k;
        end
        P_k = N0_cheng * diag(p_cheng);
        
        u_cheng = u_cheng + (theta_cheng / eta_cheng) * P_k * L_base * E_cheng;
        E_prev_cheng = E_cheng;
        
    end
end
fprintf('Simulation complete!\n');

% Averaging Errors
err_out_prop = err_out_prop / Ks; err_in_prop = err_in_prop / Ks;
err_out_shen = err_out_shen / Ks; err_in_shen = err_in_shen / Ks;
err_out_cheng = err_out_cheng / Ks; err_in_cheng = err_in_cheng / Ks;

%% 6. Plotting Results
time_axis = 0 : (N-1); 
iter_axis = 0 : (num_iter-1); 

yd_y1 = yd(1:2:end); yd_y2 = yd(2:2:end);
Y5_y1 = Y_prop_k5(1:2:end); Y5_y2 = Y_prop_k5(2:2:end);

% =========================================================================
% Figure 1: Tracking Performance
% =========================================================================
figure(1);
hold on; box on;

h1 = plot(time_axis, yd_y1, '--d', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'MarkerSize', 6);
     plot(time_axis, yd_y2, '--d', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, 'MarkerSize', 6);

h2 = plot(time_axis, Y5_y1, '-.+', 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'MarkerSize', 6);
     plot(time_axis, Y5_y2, '-.+', 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'MarkerSize', 6);

xlim([0, N-1]); 
ylim([-2, 2]); 

xlabel('Time Step \it t', 'FontSize',15);
ylabel('Trajectory Tracking','FontSize',15,'Interpreter','latex');
legend([h1, h2], {'$y_{d}(t)$','$y_{5}(t)$'}, ...
      'Location', 'southwest', 'Interpreter', 'latex');

% =========================================================================
% Figure 2: Output Error Evolution 
% =========================================================================
figure(2);
semilogy(iter_axis, err_out_shen, ':', 'Color', [0.0000, 0.4470, 0.7410]); hold on; grid on; box on;
semilogy(iter_axis, err_out_cheng, '--', 'Color', [0.9290, 0.6940, 0.1250]);
semilogy(iter_axis, err_out_prop, '-', 'Color', [0.8500, 0.3250, 0.0980]);

% xi = err_out_prop(1); 
xi = 4.2; 

semilogy(iter_axis, xi./(iter_axis+0.001).^0.72,':','color',[0 0 0]);
xlim([0, num_iter-1]); 
xlabel('Iteration Number \it k','FontSize',15);
ylabel('Tracking Error','FontSize',15,'Interpreter','latex');
legend('SA-Shen [33]', 'SA-Cheng [10]', 'Algorithm 3.1',...
    '$$O(1/k^{0.72})$$', 'Location', 'northeast', 'Interpreter', 'latex');

% =========================================================================
% Figure 3: Input Error Evolution
% =========================================================================
figure(3);
semilogy(iter_axis, err_in_shen, ':', 'Color', [0.0000, 0.4470, 0.7410]); hold on; grid on; box on;
semilogy(iter_axis, err_in_cheng, '--', 'Color', [0.9290, 0.6940, 0.1250]);
semilogy(iter_axis, err_in_prop,  '-', 'Color', [0.8500, 0.3250, 0.0980]);

xlim([0, num_iter-1]); 
xlabel('Iteration Number \it k','FontSize',15);
ylabel('Input Error','FontSize',15);
legend('SA-Shen [33]', 'SA-Cheng [10]', 'Proposed', 'Location', 'northeast');