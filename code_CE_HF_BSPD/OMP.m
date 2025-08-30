function [X] = OMP(Y, Phi, L, N, M)
%%% initialization
R = Y;
%%% OMP
X = zeros(N, M);
for m = 1:M
    support = [];
    for i = 1 : L
        t = Phi'*R(:, m);
        t = abs(t).^2;
        [~,order] = max(t);
        support = [support order];
        Phi_s = Phi(:,support);
        x = zeros(N,1);
        x(support) = pinv(Phi_s'*Phi_s)*Phi_s'*Y(:, m);
        R(:, m) = Y(:, m) - Phi*x; 
    end
    
    Phi_s = Phi(:, support);
    X(support, m) = pinv(Phi_s'*Phi_s)*Phi_s'*Y(:, m);
end