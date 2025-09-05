function [h, hc, hf, hn, r, r0, theta_f, theta_n, G] = generate_hybrid_field_channel(Nt, Ln, Lf, d, fc, B, M, Rmin, Rmax)
% Lf: far-field path num
% Ln: near-field path num
% Nt: antenna num

sector = 2*pi/3;
c = 3e8;

%far-field channel
rmax=180;
rmin=150;
r = (rand(1, Lf) * (rmax - rmin) + rmin);
theta_f = rand(1, Lf) * sector - sector/2;
gf = (rand(Lf, 1) + 1j*rand(Lf, 1))/sqrt(2);
Gf = zeros(Lf, M);
Hf = zeros( Nt, M+1 );

if Lf > 0
    for m = 1:M+1
       if m == M+1
            f = fc;
       else
            f = fc+B/M*(m-1-(M-1)/2);
       end
       lambda = c/f;
       for l = 1:Lf
           af = far_field_manifold(Nt,theta_f(l),d,lambda);
           % NLOS path
           g = gf(l) * c / f / 4 / pi / r(l);
           Hf(:, m) = Hf(:, m) + g * exp(-1j*2*pi*f*r(l)/c) * af.';
           if m <= M
               Gf(l, m) = g * exp(-1j*2*pi*f*r(l)/c); 
           end
       end 
    end
end

%near-field channel
theta_n = rand(1, Ln) * sector - sector/2;
r0 = (rand(1, Ln) * (Rmax - Rmin) + Rmin); % near-field distance
rsn = (rand(1, Ln) * (2 - 0.5) + 0.5); % near-field scattter distance
gn = (rand(Ln, 1) + 1j*rand(Ln, 1))/sqrt(2); % degradation gain
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
           %LoS path
           if l == Ln
               g =  c / f / 4 / pi / r0(l);
               [an, ~] = near_field_manifold( Nt, d, f, r0(l), theta_n(l) );
               Hn(:, m) = Hn(:, m) + g * exp(-1j*2*pi*f*r0(l)/c) * an.';
               if m <= M
                   Gn(l, m) = g * exp(-1j*2*pi*f*r0(l)/c);
               end
            % reflect path
           elseif l == Ln-1
               r_ref1 = rsn(l);
               r_ref2 = r0(l);
               r_ref = r_ref1 + r_ref2;
               g =  gn(l) * c / f / 4 / pi / r_ref;
               [an, ~] = near_field_manifold( Nt, d, f, r_ref, theta_n(l) );
               Hn(:, m) = Hn(:, m) + g * exp(-1j*2*pi*f*r_ref/c) * an.';
               if m <= M
                   Gn(l, m) = g * exp(-1j*2*pi*f*r_ref/c);
               end
               % scatter NLOS ptah
            else
               g =  gn(l) * (c / f)^2 / (4 * pi)^2 / r0(l) / rsn(l);
               [an, ~] = near_field_manifold( Nt, d, f, r0(l), theta_n(l) );
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
%Cf=isequal(hf,h);
%Cn=isequal(hn,h);
end
