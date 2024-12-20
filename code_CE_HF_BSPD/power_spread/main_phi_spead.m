clc, clear, close all

Niter = 10;
N = 512;
M=5;
fc = 300e9;
fs = 30e9;
c = 3e8;
lambda = c/fc;
d = lambda/2;

Ns = 1000;
Np = 1000;
Nq = 1000;
phi_list = linspace(-1, 1, Ns);
theta_list = linspace(-1, 1, Ns);
r0 = 3;
r = 150;
C1 = 0;
S1 = 0;
C2 = 0;
S2 = 0;
phi0 = 0.5; 
theta1 = -0.5;
fm1=285e9;
fm2=315e9;
f_mc = zeros(1,Ns);
f_mc1 = zeros(1,Ns);
gf=lambda/4/pi/r;
gn=lambda/4/pi/r0;
for iter = 1:Niter
    fprintf('iteration:[%d/%d]',iter,Niter);
    for m = 1:M
        m
        fm=fc+fs/(M)*(m-1-(M-1)/2);
        f = theta_spread(gf,fm,Ns,fc,theta_list,r0,d,N,phi0);
        f1 = phi_spread(gn,fm,Ns,fc,theta1,r0,d,N,phi_list);
        f_mc = f_mc + f;
        f_mc1 = f_mc1 + f1;
    end
f = theta_spread(gf,fm,Ns,fc,theta_list,r0,d,N,phi0);
f1 = phi_spread(gn,300e9,Ns,fc,theta0,r0,d,N,phi_list);
%f2 = phi_spread(285e9,Ns,fc,theta0,r0,d,N,phi_list);
%f3 = phi_spread(315e9,Ns,fc,theta0,r0,d,N,phi_list);
end
f1_dB = 10 * log10(abs(f1));
f_dB = 10 * log10(abs(f));
%f3_dB = 10 * log10(abs(f3));
f_mc_dB = 10 * log10(abs(f_mc)); 
f_mc1_dB = 10 * log10(abs(f_mc1));
figure;
hold on;
box on;
grid on;
plot(phi_list, f1_dB,'b--')
plot(phi_list, f_dB, 'r--')
%plot(phi_list, f3_dB, 'm-')
plot(phi_list, f_mc_dB,'r-')
plot(phi_list, f_mc1_dB,'b-')
y = ylim; 
line([0.5 0.5], y, 'Color', 'k', 'LineStyle', '--'); 
line([-0.5 -0.5], y, 'Color', 'k', 'LineStyle', '--'); 
xlabel('far-field angle $\phi_{l,m_{2}}$', 'interpreter', 'latex')
ylabel('Normalized interference power (dB)', 'interpreter', 'latex')
legend({'generated by far-field component, r=150, $\phi_{l}$=0.5, M=1', 'generated by near-field component, r=3, $\theta_{l}$=-0.5, M=1','generated by far-field component, r=150, $\phi_{l}$=0.5, M=5', 'generated by near-field component, r=3, $\theta_{l}$=-0.5, M=5'}, 'interpreter', 'latex', 'fontsize', 10);
colormap('jet')