function [Hsf_hat] = MoOMP(Z, Phi, angle_codebook, L, M, sin_theta_pattern)
    % stage 1 
    [Nt, D] = size(angle_codebook);
    %%% initialization
    R = Z;
    support = zeros(L, M);
    angle_s = zeros(Nt,L, M);
    %%% P-SOMP
    for path = 1 : L
        
        t = (Phi*angle_codebook)'*R ;
        
        pattern = zeros(D, 1);
        for m = 1:M
            sin_theta_pattern_m = sin_theta_pattern(:, m);
            pattern = pattern + abs(t ( sin_theta_pattern_m, m ) ).^2/M;
        end
%         t = sum(abs(t).^2, 3);
        
        
        [~, pattern_max] = max(pattern);
        % update
        X_s = zeros(path, M);
        R = zeros(size(Z));
        for m = 1:M
            support_m = sin_theta_pattern( pattern_max, m );
            support( path, m) = support_m;
            
            angle_s_m = angle_codebook(:, sin_theta_pattern( pattern_max, m ) );
            angle_s(:, path, m) = angle_s_m;
            
            Phi_s_m = Phi * angle_s(:, 1:path, m);
            X_s(:, m) = pinv(Phi_s_m'*Phi_s_m)*Phi_s_m'*Z(:, m);
            R(:, m) = Z(:, m) - Phi_s_m * X_s(:, m); 
        end   
    end
    Hsf_hat = zeros(Nt, M);
    for m = 1:M
        Hsf_hat(:, m) = angle_s(:, :, m) * X_s(:, m);
    end   
end

