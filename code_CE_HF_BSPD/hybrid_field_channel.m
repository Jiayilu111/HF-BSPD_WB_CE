function [h, hc, hf, hn, r0, theta_n, G] = hybrid_field_channel(Nt, Ln, Lf, d, fc, B, M, Rmin, Rmax)
% Lf: far-field path num
% Ln: near-field path num
% Nt: antenna num

sector = 2*pi/3;
c = 3e8;
rmax=180;
rmin=150;
r = (rand(1, Lf) * (rmax - rmin) + rmin);
theta_f = rand(1, Lf) * sector - sector/2;
gf = (rand(Lf, 1) + 1j*rand(Lf, 1))/sqrt(2);
Gf = zeros(Lf, M);
Hf = zeros( Nt, M+1 );
if Lf>0
    for m = 1:M+1
       if m == M+1
            f = fc;
       else
            f=fc+B/M*(m-1-(M-1)/2);
       end
       lambda = c/f;
       for l = 1:Lf
           af = far_field_manifold(Nt,theta_f(l), d, lambda);
           g = gf(l);
           Hf(:, m) = Hf(:, m) + g * exp(-1j*2*pi*f*r(l)/c) * af.';
           if m <= M
                Gf(l, m) = g ;
           end
       end 
    end
end
theta_n = rand(1, Ln) * sector - sector/2;
r0 = (rand(1, Ln) * (Rmax - Rmin) + Rmin);
nn = -(Nt-1)/2:1:(Nt-1)/2;
gn = (rand(Ln, 1) + 1j*rand(Ln, 1))/sqrt(2);
Gn = zeros(Ln, M);
Hn = zeros( Nt, M+1 );

c = 3e8;
if Ln>0
    for m = 1:M+1
       if m == M+1
            f = fc;
       else
            f=fc+B/(M)*(m-1-(M-1)/2);
       end

       for l = 1:Ln
           an = near_field_manifold( Nt, d, f, r0(l), theta_n(l) );
           g =  gn(l);
           Hn(:, m) = Hn(:, m) + g * exp(-1j*2*pi*f*r0(l)/c) * an.';
           if m <= M
                Gn(l, m) = g * exp(-1j*2*pi*f*r0(l)/c);
           end
       end 
    end
end
hcf = Hf(:,M+1);
hf = Hf(:,1:M);
hcn = Hn(:,M+1);
hn = Hn(:,1:M);
h = hf+hn;
hc = hcf+hcn;
G = [Gn;Gf];
h = h*sqrt(Nt/(Ln+Lf));
end
