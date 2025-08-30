clc,clear,close all
N_iter = 100; 
%%% system parameters
Nt= 512; % number of beams (transmit antennas)
N_RF = 4; % number of RF chains
M = 128; % number of subcarriers
L = 6; % number of paths per user
Ln = 3;
Lf = 3;
Q_list = linspace(16, 64, 9);
Q_len = length(Q_list);
fc = 300e9; % carrier frequency
Rmin = 3;
Rmax = 10;
tmax = 20e-9; % maximum of the path delay
fs = 0.1 * fc;
f = zeros(1,M);
for m = 1:M
    f(m)=fc+fs/(M)*(m-1-(M-1)/2);
end
SNR_dB =10;
SNR_linear = 10.^(SNR_dB/10.);
sigma2 = 1/SNR_linear;
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;
RD = 2 * (Nt * d)^2  / lambda_c;
eps = 1e-3;

NMSE0 = zeros(1,length(Q_list));
NMSE1 = zeros(1,length(Q_list));
NMSE2 = zeros(1,length(Q_list));
NMSE3 = zeros(1,length(Q_list));
NMSE4 = zeros(1,length(Q_list));
NMSE5 = zeros(1,length(Q_list));
NMSE6 = zeros(1,length(Q_list));
NMSE7 = zeros(1,length(Q_list));
NMSE8 = zeros(1,length(Q_list));

%DFT      
D = 2 * Nt; %字典规模
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
%Monte Carlo(time consuming)
Rh=zeros(Nt,Nt);
for s=1:10000
    [h, hc, hf, hn, r, theta, G]=hybrid_field_channel(Nt, Ln, Lf, d, fc, fs, M, Rmin, Rmax);
    Rh=Rh+h*h';
end
Rh=Rh./(10000);
parfor i_q = 1:Q_len
    Q = Q_list(i_q);
    Hfn = zeros(Nt,M);
    Hnn = zeros(Nt,M);
    % generate bilinear code
    [alpha_pattern, sin_theta_pattern] = generate_BiLinear(f, fc, alpha_list, sin_theta_list);
    
    error0 = 0; error1 = 0; error2 = 0; error3 = 0; error4 = 0; error5 = 0; error6 = 0;error7 = 0; error8 = 0;
    for iter = 1:N_iter
        fprintf('Q = %.4f [%d/%d] | iteration:[%d/%d] | run %.4f s\n', Q, i_q, Q_len, iter, N_iter, etime(clock, t0)); 
        
        % Wideband spatial channel
        [H, hc, hf, hn, ~, ~, ~, ~, G] = generate_hybrid_field_channel(Nt, Ln, Lf, d, fc, fs, M, Rmin, Rmax);
        Hsf = channel_norm(H);
        if Lf > 0
            Hfn = channel_norm(hf);
        end
        if Ln > 0
            Hnn = channel_norm(hn);
        end
        Phi = (2*randi(2,Q*N_RF,Nt) - 3)/sqrt(Nt);
        noise = sqrt(sigma2)*(randn(Q*N_RF,M)+1i*randn(Q*N_RF,M))/sqrt(2); 

        % adaptive selecting matrix
        Z = Phi*Hsf + noise;
        Znorm = norm(Z, 'fro');

        %% LS
         Haf = Hsf + ( sqrt(sigma2)*(randn(Nt,M)+1i*randn(Nt,M))/sqrt(2) );
        
        %% the MMSE
        Hsf_hat3=Rh*inv(Rh+(sigma2*eye(Nt)))*Haf; 
        
         %% Far field BSPD-OMP
        [Hsf_hat4] = MoOMP(Z, Phi, angle_codebook, 12*L, M, sin_theta_pattern);

        %% Near field BPD-OMP
        [Hsf_hat5] = BiOMP(Z, Phi, polar_codebook, 12*L, M, alpha_pattern, sin_theta_pattern);
        
        %% Hybrid field BPD-OMP
        [Hsf_hat6,Hf,Hn] = BiOMP_h(Z, Phi, polar_codebook, angle_codebook, 12*Lf, 12*Ln, M, alpha_pattern, sin_theta_pattern);
        
        %% Far field OMP
        %Hsf_hat7 = OMP(Z, Phi * DFT, 12*L, D, M);
                 
        %% Near field OMP
        Hsf_hat8 = P_OMP(Z, Phi, polar_codebook, 12*L, M);
        
        %Hsf_hat6-Hf;
        %A = isequal(Hf, Hsf_hat6);
        %B=isequal(Hfn,Hsf);
        %E = isequal(Hn, Hsf_hat6);
        %disp(E)
        %F=isequal(Hnn,Hsf);
        %disp(F)

        error0 = error0 + norm(Hfn - Hf,'fro')^2/norm(Hfn,'fro')^2;
        error1 = error1 + norm(Hnn - Hn,'fro')^2/norm(Hnn,'fro')^2;
        error2 = error2 + norm(Hsf - Haf,'fro')^2/norm(Hsf,'fro')^2;
        error3 = error3 + norm(Hsf - Hsf_hat3,'fro')^2/norm(Hsf,'fro')^2;
        error4 = error4 + norm(Hsf - Hsf_hat4,'fro')^2/norm(Hsf,'fro')^2;
        error5 = error5 + norm(Hsf - Hsf_hat5,'fro')^2/norm(Hsf,'fro')^2;
        error6 = error6 + norm(Hsf - Hsf_hat6,'fro')^2/norm(Hsf,'fro')^2;
     %   error7 = error7 + norm(Hsf - DFT*Hsf_hat7,'fro')^2/norm(Hsf,'fro')^2;
        error8 = error8 + norm(Hsf - Hsf_hat8,'fro')^2/norm(Hsf,'fro')^2;
    end
    NMSE0(i_q) = error0/N_iter;
    NMSE1(i_q) = error1/N_iter;
    NMSE2(i_q) = error2/N_iter;
    NMSE3(i_q) = error3/N_iter;
    NMSE4(i_q) = error4/N_iter;
    NMSE5(i_q) = error5/N_iter;
    NMSE6(i_q) = error6/N_iter;
  %  NMSE7(i_q) = error7/N_iter;
    NMSE8(i_q) = error8/N_iter;
end


x = 1:1:size(Q_list,2);
figure; hold on; box on; grid on;
%plot(Q_list(x),10*log10(NMSE0(x)),'->','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
%plot(Q_list(x),10*log10(NMSE1(x)),'-<','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
%plot(Q_dB_list(x),10*log10(NMSE2(x)),'k--','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
plot(Q_list(x),10*log10(NMSE3(x)),'k--','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
plot(Q_list(x),10*log10(NMSE4(x)),'-*','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
plot(Q_list(x),10*log10(NMSE5(x)),'-o','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
plot(Q_list(x),10*log10(NMSE6(x)),'r-d','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
%plot(Q_list(x),10*log10(NMSE7(x)),'-*','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');
plot(Q_list(x),10*log10(NMSE8(x)),'-<','Linewidth',1.4,'markersize',5,'MarkerFaceColor','w');


xlabel('pilot length', 'interpreter', 'latex', 'fontsize', 12);
ylabel('NMSE [dB]', 'interpreter', 'latex', 'fontsize', 12);
%legend({'LS','MMSE','FF-BPD','NF-BPD','HF-BPD','OMP','P-OMP'}, 'interpreter', 'latex', 'fontsize', 10);
legend({'MMSE','FF-BPD','NF-BPD','Proposed HF-BSPD','P-OMP'}, 'interpreter', 'latex', 'fontsize', 10);
%title('near-first');


