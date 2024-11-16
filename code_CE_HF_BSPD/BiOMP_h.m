function [Hsf_hat,Hf,Hn] = BiOMP_h(Z, Phi, polar_codebook, angle_codebook, Lf, Ln, M, alpha_pattern, sin_theta_pattern)    
    % stage 1 
    R = Z;
    [Nt, D, S] = size(polar_codebook);
    %%% initialization
    %R_n = Z;
    support_n = zeros(2, Ln, M);
    polar_s = zeros(Nt,Ln, M);

    %%% OMP_n
    for path1 = 1 : Ln
        t = zeros(D, S, M);
        for s = 1:S
            t(:, s, :) = (Phi*polar_codebook(:, :, s))'*R ;%Relevance matrix
        end
        pattern = zeros(D, S);
        for m = 1:M
            sin_theta_pattern_m = sin_theta_pattern(:, m);
            alpha_pattern_m = alpha_pattern(:, m);
            pattern = pattern + abs(t(sin_theta_pattern_m, alpha_pattern_m, m)).^2/M;
        end
        pattern_max = max(max(pattern));%(angle offset+distance offset)max relevance
        [row, col] = find( pattern_max == pattern );%corresponding to samples
        % update
        H_p = zeros(path1, M);
        R = zeros(size(Z));
        for m = 1:M
            support_m = [sin_theta_pattern(row, m), alpha_pattern(col, m)];%Minimum error support set
            support_n(:, path1, m) = support_m;
            polar_s_m = polar_codebook(:, sin_theta_pattern(row, m), alpha_pattern(col, m));
            polar_s(:, path1, m) = polar_s_m;
            Phi_n_m = Phi * polar_s(:, 1:path1, m);
            H_p(:, m) = pinv(Phi_n_m'*Phi_n_m)*Phi_n_m'*Z(:, m);
            R(:, m) = Z(:, m) - Phi_n_m * H_p(:, m); 
        end
    end
    % stage 2
    Rf=R;
    [Nt, D] = size(angle_codebook);
    %%% initialization
    support_f = zeros(Lf, M);
    angle_s = zeros(Nt, Lf, M);
    
    for path2 = 1 :Lf
        t = (Phi * angle_codebook)' * Rf;
        pattern_f = zeros(D, 1);
        for m = 1:M
            sin_theta_pattern_m = sin_theta_pattern(:, m);
            pattern_f = pattern_f + abs(t(sin_theta_pattern_m, m)).^2/M;
        end
        [~,pattern_f_max] = max(pattern_f);
        H_a = zeros(path2, M);
        for m = 1:M
            support_f_m = sin_theta_pattern(pattern_f_max, m);
            support_f(path2, m) = support_f_m;
            angle_s_m = angle_codebook(:,sin_theta_pattern(pattern_f_max, m));
            angle_s(:, path2, m) = angle_s_m;
            Phi_f_m = Phi * angle_s(:, 1:path2, m);
            H_a(:, m) = pinv(Phi_f_m' * Phi_f_m) * Phi_f_m' * R(:, m);
            Rf(:,m) = R(:, m) - Phi_f_m * H_a(:,m);
        end
    end
    Hsf_hat = zeros(Nt, M);
    Hf = zeros(Nt,M);
    Hn = zeros(Nt,M);

    if Ln > 0
        if Lf > 0
            for m = 1:M
                Hn(:,m) = polar_s(:, :, m) * H_p(:, m);
                Hf(:,m) = angle_s(:, :, m) * H_a(:, m);
                Hsf_hat(:, m) = polar_s(:, :, m) * H_p(:, m) + angle_s(:, :, m) * H_a(:, m);
            end
        else
            for m = 1:M
                Hn(:,m) = polar_s(:, :, m) * H_p(:, m);
                Hsf_hat(:, m) = polar_s(:, :, m) * H_p(:, m);
            end
        end
    else
        for m = 1:M
            Hf(:,m) = angle_s(:, :, m) * H_a(:, m);
            Hsf_hat(:, m) = angle_s(:, :, m) * H_a(:, m);
        end
    end
%    isequal(Hsf_hat,Hn);
end
