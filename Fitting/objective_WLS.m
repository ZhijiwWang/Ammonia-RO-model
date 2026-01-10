function J = objective_WLS(z, pH_grid, pH_exp, Rej_exp, Rej_err, A_exp, A_err, lambda_A, decodeParams)
% Decode z -> (p1a,p1b,p1c,p1d,f)
[p1a,p1b,p1c,p1d,f] = decodeParams(z);

% Run the physical model on pH_grid
[Rej_grid, A_grid, ~] = run_model_with(p1a,p1b,p1c,p1d,f, pH_grid);

% Interpolate model to experimental feed-pH points
Rej_fit = interp1(pH_grid, Rej_grid, pH_exp, 'linear','extrap');
A_fit   = interp1(pH_grid, A_grid,   pH_exp, 'linear','extrap');

% Weighted residuals
wR = 1./(Rej_err.^2);
wA = 1./(A_err.^2);

resR = sqrt(wR(:)) .* (Rej_fit(:) - Rej_exp(:));
resA = sqrt(wA(:)) .* (A_fit(:)   - A_exp(:));

% OPTIONAL: include permeate pH (need pH_perm data & errors)
% resPH = sqrt(wPH(:)) .* (pHperm_fit(:) - pHperm_exp(:));

% Weighted SSE (balance A term with lambda_A)
J = sum(resR.^2) + lambda_A * sum(resA.^2);

% Robustness: penalize NaN / Inf
if ~isfinite(J); J = 1e12; end
end
