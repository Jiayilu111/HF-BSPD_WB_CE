clc, clear, close all

N = 512;
M = 256;
fc = 100e9;
c = 3e8;
fs = 10e9;
lambda = c/fc;
d = lambda/2;

Ns = 1000;
Np = 1000;
Nq = 1000;
phi = 0;
r = 3;
C1 = 0;
S1 = 0;
C2 = 0;
S2 = 0;
%f = zeros(1,Ns);
for m = 1:M
    fm_list(m)=fc+fs/(M)*(m-1-(M-1)/2)
end
for m = 1:M
    theta = 0.3;
    fm=fm_list(m);
    yta = fm/fc;

    b1 = (theta - phi) * sqrt(r / d / (1 - theta^2) / yta);
    b2 = N / 2 * sqrt(yta * d * (1 - theta^2) / r);
    sum = b1 + b2;
    cha = b1 - b2;
    x1_list = linspace(0, sum, Np);
    x2_list = linspace(0, cha, Nq);
    A1=sqrt(yta*d*(1-theta^2)/2/r);
    A2=( 2 * r * yta * (theta - phi) + yta * (N - 1) * d * (1 - theta^2)) / 4 / r / A1;
    C = exp(-1i * pi * A2^2 + (1i * pi * yta * (N-1) * theta / 2) + (1i * pi * yta * (N - 1)^2 * d * (1 - theta^2)));

    C1 = integral(@(t) cos(pi/2 * t.^2), 0, sum);
    S1 = integral(@(t) sin(pi/2 * t.^2), 0, sum);
    C2 = integral(@(t) cos(pi/2 * t.^2), 0, cha);
    S2 = integral(@(t) sin(pi/2 * t.^2), 0, cha);

    C_hat = C1 - C2;
    S_hat = S1 - S2;

    f1(m) = C^N * (C_hat + 1i * S_hat) / 2 / b2;
    %f(i) = (C_hat + 1i * S_hat) / 2 / cha
end

for m = 1:M
    theta = 0;
    yta = fm_list(m)/fc;

    b1 = (theta - phi) * sqrt(r / d / (1 - theta^2) / yta);
    b2 = N / 2 * sqrt(yta * d * (1 - theta^2) / r);
    sum = b1 + b2;
    cha = b1 - b2;
    x1_list = linspace(0, sum, Np);
    x2_list = linspace(0, cha, Nq);
    A1=sqrt(yta*d*(1-theta^2)/2/r);
    A2=( 2 * r * yta * (theta - phi) + yta * (N - 1) * d * (1 - theta^2)) / 4 / r / A1;
    C = exp(-1i * pi * A2^2 + (1i * pi * yta * (N-1) * theta / 2) + (1i * pi * yta * (N - 1)^2 * d * (1 - theta^2)));

    C1 = integral(@(t) cos(pi/2 * t.^2), 0, sum);
    S1 = integral(@(t) sin(pi/2 * t.^2), 0, sum);
    C2 = integral(@(t) cos(pi/2 * t.^2), 0, cha);
    S2 = integral(@(t) sin(pi/2 * t.^2), 0, cha);

    C_hat = C1 - C2;
    S_hat = S1 - S2;

    f2(m) = C^N * (C_hat + 1i * S_hat) / 2 / b2;
    %f(i) = (C_hat + 1i * S_hat) / 2 / cha
end

f1_dB = 10 * log10(abs(f1));  
f2_dB = 10 * log10(abs(f2)); 

figure;
hold on;
box on;
grid on;
plot(fm_list, f1_dB,'b-','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w')
plot(fm_list, f2_dB, 'r-','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w')
ylim([-30, 0]);
xlabel('f(Hz)', 'interpreter', 'latex')
ylabel('Normalized interference power(dB)', 'interpreter', 'latex')
legend({'$\phi$=0, $\theta$=0.3, r=3', '$\phi$=0, $\theta$=0, r=3'}, 'interpreter', 'latex', 'fontsize', 10);
colormap('jet')