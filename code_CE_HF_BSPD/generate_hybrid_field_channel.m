function [h, hc, hf, hn, r0, theta_n, G] = generate_hybrid_field_channel(Nt, Ln, Lf, d, fc, B, M, Rmin, Rmax)
% Lf: far-field path num
% Ln: near-field path num
% Nt: antenna num

sector = 2*pi/3;
c = 3e8;
if Ln > 0 && Lf > 0
    Lnl = randi(2) - 1;
    Lfl = 1 - Lnl;
else
    Lnl = Ln > 0;
    Lfl = Lf > 0;
end

%far-field channel
rmax=180;
rmin=150;
r = (rand(1, Lf) * (rmax - rmin) + rmin);
rsf = (rand(1, Lf) * (2 - 0.5) + 0.5);
theta_f = rand(1, Lf) * sector - sector/2;
gf = (randn(Lf, 1) + 1j*randn(Lf, 1))/sqrt(2);
Gf = zeros(Lf, M);
Hf = zeros( Nt, M+1 );
if Lf > 0
    for m = 1:M+1
       if m == M+1
           f = fc;
       else
           f = fc+B/M*(m-1-(M-1)/2);
       end
       for l = 1:Lf
           af = far_field_manifold(Nt,theta_f(l));
           %LOS
           if l == Lf && Lfl == 1
               g = c / f / 4 / pi / r(l);
               Hf(:, m) = Hf(:, m) + g * exp(-1j*2*pi*f*r(l)/c) * af.';
               %Hf(:, m) = Hf(:, m) + g * af.';
               if m <= M
                  Gf(l, m) = g * exp(-1j*2*pi*f*r(l)/c); 
               end
           %NLOS
           else
               g = gf(l) * (c / f)^2 / (4 * pi)^2 / r(l) / rsf(l);
               Hf(:, m) = Hf(:, m) + g * exp(-1j*2*pi*f*(r(l)+rsf(l))/c) * af.';
               %Hf(:, m) = Hf(:, m) + g * af.';
               if m <= M
                   Gf(l, m) = g * exp(-1j*2*pi*f*(r(l)+rsf(l))/c) ;
               end
           end 
       end
    end
end

%near-field channel
theta_n = rand(1, Ln) * sector - sector/2;
r0 = (rand(1, Ln) * (Rmax - Rmin) + Rmin);
rsn = (rand(1, Ln) * (2 - 0.5) + 0.5);
nn = -(Nt-1)/2:1:(Nt-1)/2;
gn = (randn(Ln, 1) + 1j*randn(Ln, 1))/sqrt(2);
Gn = zeros(Ln, M);
Hn = zeros( Nt, M+1 );

if Ln > 0
    for m = 1:M+1
       if m == M+1
           f = fc;
       else
           f=fc+B/(M)*(m-1-(M-1)/2);
       end
       for l = 1:Ln
           an = near_field_manifold( Nt, d, f, r0(l), theta_n(l) );
           %LOS
           if l == Ln && Lnl == 1
               g =  c / f / 4 / pi / r0(l);
               Hn(:, m) = Hn(:, m) + g * exp(-1j*2*pi*f*r0(l)/c) * an.';
               if m <= M
                   Gn(l, m) = g * exp(-1j*2*pi*f*r0(l)/c);
               end
           %NLOS
           else
               g =  gn(l) * (c / f)^2 / (4 * pi)^2 / r0(l) / rsn(l);
               Hn(:, m) = Hn(:, m) + g * exp(-1j*2*pi*f*(r0(l)+rsn(l))/c) * an.';
               if m <= M
                   Gn(l, m) = g * exp(-1j*2*pi*f*(r0(l)+rsn(l))/c);
               end
           end 
        end
    end
end

%combining
hcf = Hf(:,M+1);
hf = Hf(:,1:M);
hcn = Hn(:,M+1);
hn = Hn(:,1:M);
h = hf + hn;
hc = hcf + hcn;
G = [Gn;Gf];
h = h * sqrt(Nt/(Ln+Lf));

end
