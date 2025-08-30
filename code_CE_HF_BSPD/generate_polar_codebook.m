function [w, sin_theta_list, r_list, alpha_list] = generate_polar_codebook(Nt, d, lambda_c, D, beta, rho_min, rho_max)
    c = 3e8;
    fc = c/lambda_c;
    sin_theta_list = (-1 + 2/D : 2/D : 1)';
%     D = size(sin_theta_list, 2);
    
    Z_delta = (Nt * d)^2 / 2 / lambda_c / beta^2;
    Smax = ceil(Z_delta/rho_min) ;
    Smin = floor(Z_delta/rho_max) + 1;
    S = Smax - Smin + 2;
    r_list = zeros(S,1);
    r_list(1) = 200 * Nt^2 * d^2 / lambda_c;
    r_list(2:end) = Z_delta ./ (Smin:Smax);
    
    alpha_list = 1./2./r_list;

    %% polar domain codebook
    w = zeros(Nt, D, S);
%     polar_list = zeros(D, S);
    for s = 1:S
        for idx = 1:D
           sin_theta = sin_theta_list(idx);
           alpha = alpha_list(s);
           w(:, idx, s) = polar_domain_manifold(Nt, d, fc, asin(sin_theta), alpha).';
        end 
    end

end
