function J = obj_A_only(zf, z4_fixed, pH_grid, pH_exp, A_exp, decode_p1, f_decode)
    [p1a,p1b,p1c,p1d] = decode_p1(z4_fixed);  
    f = f_decode(zf);
    [~, A_grid] = run_model_with(p1a,p1b,p1c,p1d,f, pH_grid);
    A_fit = A_grid;                           % pH_grid = pH_exp
    
    r  = A_fit(:)-A_exp(:);
    J  = sum(r.^2);
    if ~isfinite(J), J = 1e12; end
end