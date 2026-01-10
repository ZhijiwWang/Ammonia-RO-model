%% =====================  DATA  =====================
parpool;
clear; clc;

pH_exp  = [6.245, 7.035, 7.815, 8.465, 9.69];
Rej_exp = [0.88,  0.89,  0.88,  0.83,  0.33];
A_exp   = [0.5294533739, 0.5316164418, 0.554543418, 0.5181868471, 0.5632406788];

% Measurement standard deviations (used as weights 1/sigma^2)
Rej_err = [0.01860807319,0.03191890822,0.03033006379,0.001969945578,0.05535572363];
A_err   = [0.08085508874,0.08507095201,0.01049512449,0.02670081646,0.004527795146];

% Sanity: avoid divide-by-zero
Rej_err = max(Rej_err, 1e-6);
A_err   = max(A_err,   1e-6);

%% =====================  OPTIONS  =====================
% pH grid used by the model (same as your code)
pH_grid = 5.3:0.1:11;

% Bounds we want to impose (handled via smooth transforms)
p1_lb = 0.01;  p1_ub = 0.90;     % steric partition p1*
f_lb  = 1e13;  f_ub  = 1e17;     % fluid–membrane friction

% Balance weights between Rej and A terms (1 = equal after standardization)
lambda_A = 1.0;

%% =====================  PARAM TRANSFORMS  =====================
% We optimize in unconstrained space z, then map to physical parameters:
% p1 = p1_lb + (p1_ub - p1_lb) * sigmoid(z)
% f  = 10^( log10(f_lb) + (log10(f_ub)-log10(f_lb)) * sigmoid(z_f) )

sigmoid = @(z) 1./(1+exp(-z));
decodeParams = @(z) deal( ...
    p1_lb + (p1_ub-p1_lb)*sigmoid(z(1)), ... % p1a
    p1_lb + (p1_ub-p1_lb)*sigmoid(z(2)), ... % p1b
    p1_lb + (p1_ub-p1_lb)*sigmoid(z(3)), ... % p1c
    p1_lb + (p1_ub-p1_lb)*sigmoid(z(4)), ... % p1d
    10.^( log10(f_lb) + (log10(f_ub)-log10(f_lb))*sigmoid(z(5)) ) ... % f
);

% Build an initial guess z0 from your nominal values
p1a0=0.09; p1b0=0.078; p1c0=0.056; p1d0=0.01; f0=1.75e15;
enc = @(x,lb,ub) log( (x-lb) / max(1e-12,(ub-x)) ); % inverse of sigmoid scaling
z0 = zeros(1,5);
z0(1) = enc(p1a0, p1_lb, p1_ub);
z0(2) = enc(p1b0, p1_lb, p1_ub);
z0(3) = enc(p1c0, p1_lb, p1_ub);
z0(4) = enc(p1d0, p1_lb, p1_ub);
z0(5) = enc( log10(f0), log10(f_lb), log10(f_ub) ); % do the “log-space” bounds

%% =====================  OBJECTIVE (Nelder–Mead)  =====================
% Weighted least squares of [Rej, A] residuals evaluated at pH_exp.
obj = @(z) objective_WLS(z, pH_grid, pH_exp, ...
                         Rej_exp, Rej_err, A_exp, A_err, lambda_A, decodeParams);

opts = optimset('Display','iter','MaxFunEvals',5e4,'MaxIter',5e4,'TolX',1e-8,'TolFun',1e-10);

% (Optional) Multiple starts for robustness
bestJ = inf; bestZ = z0;
for s = 1:5
    if s==1
        z_init = z0;
    else
        z_init = z0 + 0.2*randn(size(z0));
    end
    [z_hat, J_hat] = fminsearch(obj, z_init, opts);
    if J_hat < bestJ
        bestJ = J_hat; bestZ = z_hat;
    end
end

% Decode best-fit params
[p1a_hat,p1b_hat,p1c_hat,p1d_hat,f_hat] = decodeParams(bestZ);

% Evaluate model on grid and at experimental pH for reporting
[Rej_grid, A_grid, pHperm_grid] = run_model_with(p1a_hat,p1b_hat,p1c_hat,p1d_hat,f_hat, pH_grid);
Rej_fit_at_exp = interp1(pH_grid, Rej_grid, pH_exp, 'linear','extrap');
A_fit_at_exp   = interp1(pH_grid, A_grid,   pH_exp, 'linear','extrap');

%% =====================  REPORT  =====================
fprintf('\n==== Best-fit parameters (Nelder–Mead) ====\n');
fprintf('p1a = %.5f\n', p1a_hat);
fprintf('p1b = %.5f\n', p1b_hat);
fprintf('p1c = %.5f\n', p1c_hat);
fprintf('p1d = %.5f\n', p1d_hat);
fprintf('f   = %.3e (mol·s·m^-5)\n', f_hat);
fprintf('Objective (weighted SSE) = %.6g\n', bestJ);

% Simple goodness-of-fit summary
wR = 1./(Rej_err.^2);
wA = 1./(A_err.^2);
RSS_Rej = sum( wR .* (Rej_fit_at_exp - Rej_exp).^2 );
RSS_A   = sum( wA .* (A_fit_at_exp   - A_exp  ).^2 );
fprintf('RSS_Rej = %.6g,   RSS_A = %.6g\n', RSS_Rej, RSS_A);

% Quick plot
figure; 
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot(pH_grid, Rej_grid,'-','LineWidth',1.8);
errorbar(pH_exp, Rej_exp, Rej_err,'o','MarkerFaceColor',[.2 .6 .9]);
xlabel('feed pH'); ylabel('Rejection');
title('Fit on Rej'); legend('Model','Exp','Location','best'); grid on;

nexttile; hold on;
plot(pH_grid, A_grid,'-','LineWidth',1.8);
errorbar(pH_exp, A_exp, A_err,'s','MarkerFaceColor',[.9 .4 .2]);
xlabel('feed pH'); ylabel('A');
title('Fit on A'); legend('Model','Exp','Location','best'); grid on;

%% =====================  BOOTSTRAP CIs  =====================
% Assumes you already have:
% - bestZ (best params in unconstrained space)
% - decodeParams, objective_WLS, run_model_with, pH_grid, pH_exp
% - Rej_exp, Rej_err, A_exp, A_err, lambda_A, opts

rng(2025);                     % reproducible
B = 500;                       % e.g., 500 (use 1000 for final)
theta_boot = nan(B,5);         % [p1a p1b p1c p1d f] in original scales

% Fitted model values at experimental pH (baseline for residual bootstrap)
[p1a_hat,p1b_hat,p1c_hat,p1d_hat,f_hat] = decodeParams(bestZ);
[Rej_grid, A_grid] = run_model_with(p1a_hat,p1b_hat,p1c_hat,p1d_hat,f_hat, pH_grid);
Rej_fit = interp1(pH_grid, Rej_grid, pH_exp, 'linear','extrap');
A_fit   = interp1(pH_grid, A_grid,   pH_exp, 'linear','extrap');

% Raw residuals (not weighted) for resampling
res_Rej = Rej_exp(:) - Rej_fit(:);
res_A   = A_exp(:)   - A_fit(:);

% Bootstrap loop
parfor b = 1:B
    try
        % --- residual resampling (nonparametric) ---
        idxR = randi(numel(res_Rej), [numel(res_Rej),1]);
        idxA = randi(numel(res_A),   [numel(res_A),1]);
        Rej_syn = Rej_fit(:) + res_Rej(idxR);
        A_syn   = A_fit(:)   + res_A(idxA);

        % Build a temporary objective that uses synthetic data
        obj_b = @(z) objective_WLS(z, pH_grid, pH_exp, ...
                            Rej_syn, Rej_err, A_syn, A_err, lambda_A, decodeParams);

        % Fit starting from bestZ (robust & fast)
        z0b = bestZ + 0.05*randn(size(bestZ));     % slight jitter helps robustness
        [zb, Jb] = fminsearch(obj_b, z0b, opts);   

        % Save in original scales
        [p1a_b,p1b_b,p1c_b,p1d_b,f_b] = decodeParams(zb);
        theta_boot(b,:) = [p1a_b,p1b_b,p1c_b,p1d_b,f_b];

    catch
        % If a rare failure happens, leave row as NaN (we'll ignore it)
        continue;
    end
end

% Remove failed runs
theta_boot = theta_boot(all(isfinite(theta_boot),2),:);

% Percentile CIs
pcts = [2.5 50 97.5];
CI_p1a = prctile(theta_boot(:,1), pcts);
CI_p1b = prctile(theta_boot(:,2), pcts);
CI_p1c = prctile(theta_boot(:,3), pcts);
CI_p1d = prctile(theta_boot(:,4), pcts);
CI_f   = prctile(theta_boot(:,5), pcts);

fprintf('\n==== Bootstrap 95%% Confidence Intervals (percentile) ====\n');
fprintf('p1a: median=%.5f  [%.5f, %.5f]\n', CI_p1a(2), CI_p1a(1), CI_p1a(3));
fprintf('p1b: median=%.5f  [%.5f, %.5f]\n', CI_p1b(2), CI_p1b(1), CI_p1b(3));
fprintf('p1c: median=%.5f  [%.5f, %.5f]\n', CI_p1c(2), CI_p1c(1), CI_p1c(3));
fprintf('p1d: median=%.5f  [%.5f, %.5f]\n', CI_p1d(2), CI_p1d(1), CI_p1d(3));
fprintf('f  : median=%.3e  [%.3e, %.3e]\n', CI_f(2), CI_f(1), CI_f(3));

% (Optional) quick visualization of bootstrap distributions
figure; tiledlayout(1,5,'TileSpacing','compact','Padding','compact');
lbls = {'p1a','p1b','p1c','p1d','f'};
for k=1:5
    nexttile; histogram(theta_boot(:,k), 30); title(lbls{k}); grid on;
end

