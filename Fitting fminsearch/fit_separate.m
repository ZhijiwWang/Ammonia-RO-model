%% =====================  DATA  =====================

clear; clc;
if isempty(gcp('nocreate'))
    parpool('local', 4); % paralell pool
end

pH_exp  = [6.245, 7.035, 7.815, 8.465, 9.69];
Rej_exp = [0.88,  0.89,  0.88,  0.83,  0.33];
A_exp   = [0.5294533739, 0.5316164418, 0.554543418, 0.5181868471, 0.5632406788];
SP_exp = [0.7190, 0.6127, 0.7563, 1.0304, 10.5331];

% Measurement standard deviations (used as weights 1/sigma^2)
Rej_err = [0.01860807319,0.03191890822,0.03033006379,0.001969945578,0.05535572363];
A_err   = [0.08085508874,0.08507095201,0.01049512449,0.02670081646,0.004527795146];
SP_err = [0.1333335644,0.2283048028,0.2220197615,0.142412154,4.067208946];
SP_err = max(SP_err, 1e-6);   % avoid divide 0
lambda_SP = 0.05;    
% Sanity: avoid divide-by-zero
Rej_err = max(Rej_err, 1e-6);
A_err   = max(A_err,   1e-6);

% --- grids ---
pH_grid = pH_exp;   % only calulate the pH point with experimental data to accelerate computing

% --- bounds ---
p1_lb = 0.04;  p1_ub = 0.2;      % p1a–p1d in (0,0.2) 
p1_lam = 0.005; p1_uam = 0.015;
f_lb  = 1e13;   f_ub  = 1e17;

% --- transforms ---
sigmoid   = @(z) 1./(1+exp(-z));
inv_sig   = @(x) log( (x)./(1-x) );                         % for (0,1)
map01     = @(z,lb,ub) lb + (ub-lb)*sigmoid(z);
inv_map01 = @(x,lb,ub) log( (x-lb)./max(1e-12,(ub-x)) );    % inverse

f_decode = @(z) 10.^( log10(f_lb) + (log10(f_ub)-log10(f_lb))*sigmoid(z) );

decode_p1 = @(z4) deal( ...
    map01(z4(1),p1_lb,p1_ub), ...
    map01(z4(2),p1_lb,p1_ub), ...
    map01(z4(3),p1_lb,p1_ub), ...
    map01(z4(4),p1_lam,p1_uam) );

% --- initial guesses -> unconstrained space ---
p1a0=0.1; p1b0=0.09; p1c0=0.05; p1d0=0.006; f0=1.75e15;

z_p1_0 = [inv_map01(p1a0,p1_lb,p1_ub), ...
          inv_map01(p1b0,p1_lb,p1_ub), ...
          inv_map01(p1c0,p1_lb,p1_ub), ...
          inv_map01(p1d0,p1_lam,p1_uam)];

f_mid   = (log10(f0)-log10(f_lb))/(log10(f_ub)-log10(f_lb));
z_f_0   = inv_sig( max(1e-9,min(1-1e-9,f_mid)) );

opts = optimset('Display','iter', ...
    'TolFun', 1e-5, ...   
    'TolX',   1e-5, ...  
    'MaxIter', 5e2, ...
    'MaxFunEvals', 5e2);

% --------- Block 1:  p1a–p1d（f fixed）---------
nStarts = 4;
J1_all   = inf(1, nStarts);                 
Z1_all   = nan(nStarts, numel(z_p1_0));    

parfor s = 1:nStarts
    z_init = z_p1_0 + 0.2*randn(size(z_p1_0));

    obj1 = @(z4) obj_RejSP_only(z4, z_f_0, pH_grid, pH_exp, ...
                                Rej_exp, Rej_err, SP_exp, SP_err, lambda_SP, ...
                                decode_p1, f_decode);

    [z4_hat, J1] = fminsearch(obj1, z_init, opts);

    J1_all(s)   = J1;          
    Z1_all(s,:) = z4_hat(:).'; 
end

% select best from parfor 
[bestJ1, idx1] = min(J1_all);
best_z_p1      = Z1_all(idx1,:);

% --------- Block 2:fit f（fix p1a–p1d）---------
nStarts = 4;
J2_all = inf(1, nStarts);
Z2_all = nan(nStarts, numel(z_f_0));



parfor s = 1:nStarts
    z_init = z_f_0 + 0.2*randn(size(z_f_0));

    obj2 = @(zf) obj_A_only(zf, best_z_p1, pH_grid, pH_exp, ...
                            A_exp, decode_p1, f_decode);

    [zf_hat, J2] = fminsearch(obj2, z_init, opts);

    J2_all(s)   = J2;
    Z2_all(s,:) = zf_hat(:).';
end

[bestJ2, idx2] = min(J2_all);
best_z_f       = Z2_all(idx2,:);

[p1a_hat,p1b_hat,p1c_hat,p1d_hat] = decode_p1(best_z_p1);
f_hat = f_decode(best_z_f);


[Rej_fit_at_exp, A_fit_at_exp, pHperm] = run_model_with(p1a_hat,p1b_hat,p1c_hat,p1d_hat,f_hat, pH_grid);



RSS_Rej = sum( (Rej_fit_at_exp - Rej_exp).^2 );
RSS_A   = sum( (A_fit_at_exp   - A_exp  ).^2 );

fprintf('\n==== Block NM fit ====\n');
fprintf('p1a=%.5f, p1b=%.5f, p1c=%.5f, p1d=%.5f,  f=%.3e\n', ...
         p1a_hat,p1b_hat,p1c_hat,p1d_hat,f_hat);
fprintf('RSS_Rej=%.6g, RSS_A=%.6g\n', RSS_Rej, RSS_A);


