clc, clear, close all

N = 512;
f = 100e9;
c = 3e8;
lambda = c/f;
d = lambda/2;
%theta = 0.3;
Ns = 1000;
Np = 1000;
Nq = 1000;
%phi_list = linspace(-1, 1, Ns);
theta_list = linspace(-1, 1, Ns);
r = 10;
C1 = 0;
S1 = 0;
C2 = 0;
S2 = 0;
%f = zeros(1,Ns);


for i = 1:Ns
    phi=0.4;
    theta = theta_list(i);
    b1 = (theta - phi) * sqrt(r / d / (1 - theta^2));
    b2 = N / 2 * sqrt(d * (1 - theta^2) / r);
    sum = b1 + b2;
    cha = b1 - b2;
    x1_list = linspace(0, sum, Np);
    x2_list = linspace(0, cha, Nq);
    A1=sqrt(d*(1-theta^2)/2/r);
    A2=( 2 * r * (theta - phi) + (N - 1) * d * (1 - theta^2)) / 4 / r / A1;
    C = exp(-1i * pi * A2^2 + (1i * pi * (N-1) * theta / 2) + (1i * pi * (N - 1)^2 * d * (1 - theta^2)));


    C1 = integral(@(t) cos(pi/2 * t.^2), 0, sum);
    S1 = integral(@(t) sin(pi/2 * t.^2), 0, sum);

    C2 = integral(@(t) cos(pi/2 * t.^2), 0, cha);
    S2 = integral(@(t) sin(pi/2 * t.^2), 0, cha);

    C_hat = C1 - C2;
    S_hat = S1 - S2;

    f1(i) = C^N * (C_hat + 1i * S_hat) / 2 / b2
    %f(i) = (C_hat + 1i * S_hat) / 2 / cha
end

for i = 1:Ns
    phi = 0;
    theta = theta_list(i);
    b1 = (theta - phi) * sqrt(r / d / (1 - theta^2));
    b2 = N / 2 * sqrt(d * (1 - theta^2) / r);
    sum = b1 + b2;
    cha = b1 - b2;
    x1_list = linspace(0, sum, Np);
    x2_list = linspace(0, cha, Nq);
    A1=sqrt(d*(1-theta^2)/2/r);
    A2=( 2 * r * (theta - phi) + (N - 1) * d * (1 - theta^2)) / 4 / r / A1;
    C = exp(-1i * pi * A2^2 + (1i * pi * (N-1) * theta / 2) + (1i * pi * (N - 1)^2 * d * (1 - theta^2)));


    C1 = integral(@(t) cos(pi/2 * t.^2), 0, sum);
    S1 = integral(@(t) sin(pi/2 * t.^2), 0, sum);

    C2 = integral(@(t) cos(pi/2 * t.^2), 0, cha);
    S2 = integral(@(t) sin(pi/2 * t.^2), 0, cha);

    C_hat = C1 - C2;
    S_hat = S1 - S2;

    f2(i) = C^N * (C_hat + 1i * S_hat) / 2 / b2
    %f(i) = (C_hat + 1i * S_hat) / 2 / cha
end

f1_dB = 10 * log10(abs(f1));  
f2_dB = 10 * log10(abs(f2));  

figure; hold on; box on; grid on;
plot(theta_list, f1_dB, 'b-')
plot(theta_list, f2_dB, 'r-')
plot([0.4,0.4],[-60,20],'k--'); 
plot([0,0],[-60,20],'k--'); 
xlabel('$\theta$', 'interpreter', 'latex')
ylabel('f (dB)', 'interpreter', 'latex')
legend({'$\phi$=0.4, r=10', '$\phi$=0, r=10'}, 'interpreter', 'latex', 'fontsize', 10);
colormap('jet')