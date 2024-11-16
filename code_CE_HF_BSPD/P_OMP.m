function Hsf_hat = P_OMP(Z, Phi, polar_codebook, L, M)

    % stage 1 
    [Nt, D, S] = size(polar_codebook);
    polar_codebook2 = reshape(polar_codebook, [Nt, D*S]);
    [P,~] = size(Phi);
    Phi2 = Phi*polar_codebook2;

    %%% initialization
    R = Z;
    Hsf_hat = zeros(Nt, M);
    support = zeros(M, L);
    

    %%% P-SOMP
    for i = 1 : L
        t = Phi2'*R;
        t = abs(t).^2;
        for m = 1:M
            t_m = t(:, m);
            [~, order] = max(t_m);
            % update
            support(m, i) = order;
%             polar_s = [polar_s, polar_codebook(:, row, col)];
%             Phi_s = Phi * polar_s;
            Phi_s = Phi2(:, support(m, 1:i));
           % X_s  =  pinv(Phi_s'*Phi_s)*Phi_s'*Z(:, m);
            X_s  =  pinv(Phi_s)*Z(:, m);
            R(:, m) = Z(:, m) - Phi_s * X_s;    
        end
    end
    
    for m = 1:M
        Phi_s = Phi2(:, support(m, 1:L));
        polar_s = polar_codebook2(:, support(m, 1:L));
        %X_s  =  pinv(Phi_s'*Phi_s)*Phi_s'*Z(:, m);
        X_s  =  pinv(Phi_s)*Z(:, m);
        Hsf_hat(:, m) = polar_s * X_s;
    end
end

