function [f]=theta_spread(g,fm,Ns,fc,theta_list,r,d,N,phi0)
    for i = 1:Ns
        yta=fc/fm;
        phi = yta * phi0;
        theta = theta_list(i);
        b1 = (theta - phi) * sqrt(r / d / (1 - theta^2));
        b2 = N / 2 * sqrt(d * (1 - theta^2) / r);
        sum = b1 + b2;
        cha = b1 - b2;

        A1=sqrt(d*(1-theta^2)/2/r);
        A2=( 2 * r * (theta - phi) + (N - 1) * d * (1 - theta^2)) / 4 / r / A1;
        C = exp(-1i * pi * A2^2 + (1i * pi * (N-1) * theta / 2) + (1i * pi * (N - 1)^2 * d * (1 - theta^2)));


        C1 = integral(@(t) cos(pi/2 * t.^2), 0, sum);
        S1 = integral(@(t) sin(pi/2 * t.^2), 0, sum);

        C2 = integral(@(t) cos(pi/2 * t.^2), 0, cha);
        S2 = integral(@(t) sin(pi/2 * t.^2), 0, cha);

        C_hat = C1 - C2;
        S_hat = S1 - S2;

         f(i) = g * C^N * (C_hat + 1i * S_hat) / 2 / b2;
        %f(i) = (C_hat + 1i * S_hat) / 2 / cha
    end

end