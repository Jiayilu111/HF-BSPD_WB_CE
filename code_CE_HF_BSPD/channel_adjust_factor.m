clc,clear,close all
N_iter = 100;
%%% system parameters
Nt = 512; % number of beams (transmit antennas)
N_RF = 4; % number of RF chains
M = 128; % number of subcarriers
L = 6; % number of paths per user
a_list=[0 0.2 0.4 0.6 0.8 1];
Q=32;
fc = 300e9; % carrier frequency
Rmin = 3;
Rmax = 10;
tmax = 20e-9; % maximum of the path delay
fs = 0.1 * fc;
f = zeros(1,M);

for m = 1:M
    f(m)=fc+fs/(M)*(m-1-(M-1)/2);
end
SNR_dB = 10;
SNR_linear = 10.^(SNR_dB / 10);
sigma2 = 1/SNR_linear;
a_len = length(a_list);
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;
RD = 2 * (Nt * d)^2  / lambda_c;
eps = 1e-3;

NMSE1 = zeros(1,length(a_list));
NMSE2 = zeros(1,length(a_list));
NMSE3 = zeros(1,length(a_list));
NMSE4 = zeros(1,length(a_list));
NMSE5 = zeros(1,length(a_list));

%DFT      
D = 2 * Nt; 
row = (-(Nt - 1)/2:(Nt - 1)/2)' ;
col = -1 + 2/D : 2/D : 1 ;
DFT  =  exp( 1j*  pi * row * col );

beta = 2.5;
rho_min = 3;
rho_max = 64;
[polar_codebook, sin_theta_list, r_list, alpha_list] = generate_polar_codebook(Nt, d, lambda_c, D, beta, rho_min, rho_max);
angle_codebook = polar_codebook(:,:,1);
S = size(r_list, 1);

t0 = clock;

parfor i_l = 1:a_len
    a = a_list(i_l);
    Ln = ceil(a * L);
    Lf = L - Ln;
    Rh=zeros(Nt,Nt);
    for s=1:10000
        [h, hc, hf, hn, r, theta, G]=hybrid_field_channel(Nt, Ln, Lf, d, fc, fs, M, Rmin, Rmax);
        Rh=Rh+h*h';
    end
    Rh=Rh./(10000);
    % generate bilinear code
    [alpha_pattern, sin_theta_pattern] = generate_BiLinear(f, fc, alpha_list, sin_theta_list);
    
    error0 = 0; error1 = 0; error2 = 0; error3 = 0; error4 = 0;error5 = 0; error8 = 0; error7 = 0;
    for iter = 1:N_iter
        fprintf('a = %.4f [%d/%d] | iteration:[%d/%d] | run %.4f s\n', a, i_l, a_len, iter, N_iter, etime(clock, t0)); 
        
        % Wideband spatial channel
        [H, hc, hf, hn, ~, ~, ~, ~, G] = generate_hybrid_field_channel(Nt, Ln, Lf, d, fc, fs, M, Rmin, Rmax);
        Hsf = channel_norm(H);
        Phi = (2*randi(2,Q*N_RF,Nt) - 3)/sqrt(Nt);
        noise = sqrt(sigma2)*(randn(Q*N_RF,M)+1i*randn(Q*N_RF,M))/sqrt(2); 

        % adaptive selecting matrix
        Z = Phi*Hsf + noise;
        Znorm = norm(Z, 'fro');

        %% Far field BSPD-OMP
        [Hsf_hat1] = MoOMP(Z, Phi, angle_codebook,12*L, M, sin_theta_pattern);
        
        %% Near field BPD-OMP
        [Hsf_hat2] = BiOMP(Z, Phi, polar_codebook, 12*L, M, alpha_pattern, sin_theta_pattern);
        
        %% LS
         Haf = Hsf + ( sqrt(sigma2)*(randn(Nt,M)+1i*randn(Nt,M))/sqrt(2) );

         %% Hybrid field BPD-OMP
        [Hsf_hat3,Hf,Hn] = BiOMP_h(Z, Phi, polar_codebook, angle_codebook, 12*Lf, 12*Ln, M, alpha_pattern, sin_theta_pattern);
        
         %% the MMSE
        Hsf_hat4=Rh*inv(Rh+(sigma2*eye(Nt)))*Haf; 
        
        error1 = error1 + norm(Hsf - Hsf_hat1,'fro')^2/norm(Hsf,'fro')^2;
        error2 = error2 + norm(Hsf - Hsf_hat2,'fro')^2/norm(Hsf,'fro')^2;
        error3 = error3 + norm(Hsf - Hsf_hat3,'fro')^2/norm(Hsf,'fro')^2;
        error4 = error4 + norm(Hsf - Hsf_hat4,'fro')^2/norm(Hsf,'fro')^2;
        error5 = error5 + norm(Hsf - Haf,'fro')^2/norm(Hsf,'fro')^2;
    end
    
    NMSE1(i_l) = error1/N_iter;
    NMSE2(i_l) = error2/N_iter;
    NMSE3(i_l) = error3/N_iter;   
    NMSE4(i_l) = error4/N_iter;
    NMSE5(i_l) = error5/N_iter;
end


x = 1:1:size(a_list,2);
figure; hold on; box on; grid on;

plot(a_list(x),10*log10(NMSE1(x)),'b->','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
plot(a_list(x),10*log10(NMSE2(x)),'m-<','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
plot(a_list(x),10*log10(NMSE3(x)),'r-d','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
plot(a_list(x),10*log10(NMSE4(x)),'k--','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
%plot(a_list(x),10*log10(NMSE5(x)),'k--','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');


xlim([min(a_list(x)), max(a_list(x))])
xlabel('channel adjust factor', 'interpreter', 'latex', 'fontsize', 12);
ylabel('NMSE [dB]', 'interpreter', 'latex', 'fontsize', 12);
legend({'BSPD','BPD','HF-BPD','MMSE'}, 'interpreter', 'latex', 'fontsize', 10);