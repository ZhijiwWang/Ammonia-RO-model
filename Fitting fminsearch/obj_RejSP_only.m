function J = obj_RejSP_only(z4, zf_fixed, pH_grid, pH_exp, ...
                            Rej_exp, Rej_err, SP_exp, SP_err, lambda_SP, ...
                            decode_p1, f_decode)
    
    [p1a,p1b,p1c,p1d] = decode_p1(z4);
    f = f_decode(zf_fixed);

   
    [Rej_grid, ~, ~, SP_grid] = run_model_with(p1a,p1b,p1c,p1d,f, pH_grid);

    
    Rej_fit = Rej_grid(:);
    SP_fit  = SP_grid(:);

    
    rR = (Rej_fit - Rej_exp(:)) ;

    
    rS = (SP_fit  - SP_exp(:))  ;

    
    J = sum(rR.^2) + lambda_SP * sum(rS.^2);

    if ~isfinite(J), J = 1e12; end
end
