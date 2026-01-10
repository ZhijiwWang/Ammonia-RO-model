if isempty(gcp('nocreate')), parpool; end

% sa_sobol.m
clear; clc; rng(2025);

% ---------- Parameters for SA ----------
names = {'p1a','p1b','p1c','p1d',...
         'Kcna','Kccl','Kcam','Kcamm',...
         'e','Lm_nm','k',...
         'Kdna','Kdcl','Kdam','Kdamm' };

% +- 25% range：
lb = [0.0675, 0.0585, 0.042, 0.0075, ...      % p1a-d (dielectric/partition)
      0.15, 0.15, 0.975, 0.99, ...      % Kc 
      0.075, 112.5,  2.2e-5...     % e, Lm[nm], k[m/s]
      0.15, 0.15, 0.0495, 0.278];      % Kd      
ub = [0.1125, 0.0975, 0.07, 0.0125, ...
      0.25, 0.25, 2.0, 2.0, ...
      0.125, 187.5, 3.66e-5...
      0.15, 0.15, 0.0825, 0.464];

p = numel(names);

% ---------- Saltelli sampling ----------
Nbase = 5000;          
sA = sobolset(p,'Skip',1e3,'Leap',1e2); sA = scramble(sA,'MatousekAffineOwen');
sB = sobolset(p,'Skip',2e6,'Leap',1e2); sB = scramble(sB,'MatousekAffineOwen');

A01 = net(sA, Nbase);
B01 = net(sB, Nbase);

A = bsxfun(@plus, lb, bsxfun(@times, A01, (ub-lb)));
B = bsxfun(@plus, lb, bsxfun(@times, B01, (ub-lb)));

% ----------  A、B ----------
YA = zeros(Nbase,1); YB = zeros(Nbase,1);
parfor n=1:Nbase
    outA = eval_model(A(n,:));
    outB = eval_model(B(n,:));
    YA(n) = outA.R_mean;     
    YB(n) = outB.R_mean;
end

% ----------  A_Bi  ----------
S1  = zeros(p,1);
ST  = zeros(p,1);
for i=1:p
    ABi = A;
    ABi(:,i) = B(:,i); 
    YABi = zeros(Nbase,1);
    parfor n=1:Nbase
        out = eval_model(ABi(n,:));
        YABi(n) = out.R_mean;
    end

    % Saltelli 
    VY = var([YA; YB], 1);  
    S1(i) = (1/Nbase) * sum( YB .* (YABi - YA) ) / VY;
    ST(i) = (1/(2*Nbase)) * sum( (YA - YABi).^2 ) / VY;
    fprintf('%-8s  S1=%.3f  ST=%.3f\n', names{i}, S1(i), ST(i));
end

% ---------- plot ----------
figure; 
tiledlayout(1,2,"Padding","compact","TileSpacing","compact");
nexttile; bar(S1); title('Sobol S_1 (R\_mean)'); xticklabels(names); xtickangle(45); ylim([0 1]);
nexttile; bar(ST); title('Sobol S_T (R\_mean)'); xticklabels(names); xtickangle(45); ylim([0 1]);


S1_CI = prctile(S1_boot',[2.5 97.5])';
ST_CI = prctile(ST_boot',[2.5 97.5])';

% print
fprintf('\n=== Sobol indices with ~95%% bootstrap CIs ===\n');
for i=1:p
    fprintf('%-8s  S1=%.3f [% .3f, %.3f]   ST=%.3f [% .3f, %.3f]\n', ...
        names{i}, S1(i), S1_CI(i,1), S1_CI(i,2), ST(i), ST_CI(i,1), ST_CI(i,2));
end
