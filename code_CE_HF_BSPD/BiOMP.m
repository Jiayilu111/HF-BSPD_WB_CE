function [Hsf_hat] = BiOMP(Z, Phi, polar_codebook, L, M, alpha_pattern, sin_theta_pattern)
    % stage 1 
    [Nt, D, S] = size(polar_codebook);
    %%% initialization
    R = Z;
    support = zeros(2, L, M);
    polar_s = zeros(Nt,L, M);
    %%% P-SOMP
    for path = 1 : L
        t = zeros(D, S, M);
        for s = 1:S
            t(:, s, :) = (Phi*polar_codebook(:, :, s))'*R ;
        end
        
        pattern = zeros(D, S);
        for m = 1:M
            sin_theta_pattern_m = sin_theta_pattern(:, m);
            alpha_pattern_m = alpha_pattern(:, m);
            pattern = pattern + abs(t(sin_theta_pattern_m, alpha_pattern_m, m)).^2/M;
        end
        pattern_max = max(max(pattern));
        [row, col] = find( pattern_max == pattern );
        row=row(1);
        col=col(1);
        % update
        X_s = zeros(path, M);
        R = zeros(size(Z));
        for m = 1:M
            support_m = [sin_theta_pattern(row, m), alpha_pattern(col, m)];
            support(:, path, m) = support_m;
            polar_s_m = polar_codebook(:, sin_theta_pattern(row, m), alpha_pattern(col, m));
            polar_s(:, path, m) = polar_s_m;
            Phi_s_m = Phi * polar_s(:, 1:path, m);
            X_s(:, m) = pinv(Phi_s_m'*Phi_s_m)*Phi_s_m'*Z(:, m);
            R(:, m) = Z(:, m) - Phi_s_m * X_s(:, m);
        end   
    end
    Hsf_hat = zeros(Nt, M);
    for m = 1:M
        Hsf_hat(:, m) = polar_s(:, :, m) * X_s(:, m);
    end   
end