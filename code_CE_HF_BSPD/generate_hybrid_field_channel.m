function [h, hc, hf, hn, r0, theta_n, G] = generate_hybrid_field_channel(Nt, Ln, Lf, d, fc, B, M, Rmin, Rmax)
% Lf: far-field path num
% Ln: near-field path num
% Nt: antenna num

sector = 2*pi/3;
c = 3e8;
<<<<<<< HEAD
if Ln > 0 && Lf > 0
    Lnl = randi(2) - 1;
    Lfl = 1 - Lnl;
else
    Lnl = Ln > 0;
    Lfl = Lf > 0;
end
=======
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3

%far-field channel
rmax=180;
rmin=150;
r = (rand(1, Lf) * (rmax - rmin) + rmin);
<<<<<<< HEAD
rsf = (rand(1, Lf) * (2 - 0.5) + 0.5);
theta_f = rand(1, Lf) * sector - sector/2;
gf = (randn(Lf, 1) + 1j*randn(Lf, 1))/sqrt(2);
=======
rsf = (rand(1, Lf-1) * (2 - 0.5) + 0.5);
theta_f = rand(1, Lf) * sector - sector/2;
gf = (randn(Lf-1, 1) + 1j*randn(Lf-1, 1))/sqrt(2);
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3
Gf = zeros(Lf, M);
Hf = zeros( Nt, M+1 );
if Lf > 0
    for m = 1:M+1
       if m == M+1
<<<<<<< HEAD
           f = fc;
       else
           f = fc+B/M*(m-1-(M-1)/2);
=======
            f = fc;
       else
            f = fc+B/M*(m-1-(M-1)/2);
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3
       end
       for l = 1:Lf
           af = far_field_manifold(Nt,theta_f(l));
           %LOS
<<<<<<< HEAD
           if l == Lf && Lfl == 1
               g = c / f / 4 / pi / r(l);
               Hf(:, m) = Hf(:, m) + g * exp(-1j*2*pi*f*r(l)/c) * af.';
               %Hf(:, m) = Hf(:, m) + g * af.';
               if m <= M
                  Gf(l, m) = g * exp(-1j*2*pi*f*r(l)/c); 
=======
           if l == Lf
               g = c / f / 4 / pi / r(l);
               Hf(:, m) = Hf(:, m) + g * exp(-1j*2*pi*f*r(l)/c) * af.';
                   %Hf(:, m) = Hf(:, m) + g * af.';
               if m <= M
                   Gf(l, m) = g * exp(-1j*2*pi*f*r(l)/c); 
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3
               end
           %NLOS
           else
               g = gf(l) * (c / f)^2 / (4 * pi)^2 / r(l) / rsf(l);
               Hf(:, m) = Hf(:, m) + g * exp(-1j*2*pi*f*(r(l)+rsf(l))/c) * af.';
               %Hf(:, m) = Hf(:, m) + g * af.';
               if m <= M
<<<<<<< HEAD
                   Gf(l, m) = g * exp(-1j*2*pi*f*(r(l)+rsf(l))/c) ;
               end
           end 
       end
=======
                    Gf(l, m) = g * exp(-1j*2*pi*f*(r(l)+rsf(l))/c) ;
               end
           end
       end 
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3
    end
end

%near-field channel
theta_n = rand(1, Ln) * sector - sector/2;
r0 = (rand(1, Ln) * (Rmax - Rmin) + Rmin);
<<<<<<< HEAD
rsn = (rand(1, Ln) * (2 - 0.5) + 0.5);
nn = -(Nt-1)/2:1:(Nt-1)/2;
gn = (randn(Ln, 1) + 1j*randn(Ln, 1))/sqrt(2);
=======
rsn = (rand(1, Ln-1) * (2 - 0.5) + 0.5);
nn = -(Nt-1)/2:1:(Nt-1)/2;
gn = (randn(Ln-1, 1) + 1j*randn(Ln-1, 1))/sqrt(2);
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3
Gn = zeros(Ln, M);
Hn = zeros( Nt, M+1 );

if Ln > 0
    for m = 1:M+1
       if m == M+1
<<<<<<< HEAD
           f = fc;
       else
           f=fc+B/(M)*(m-1-(M-1)/2);
=======
            f = fc;
       else
            f=fc+B/(M)*(m-1-(M-1)/2);
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3
       end
       for l = 1:Ln
           an = near_field_manifold( Nt, d, f, r0(l), theta_n(l) );
           %LOS
<<<<<<< HEAD
           if l == Ln && Lnl == 1
               g =  c / f / 4 / pi / r0(l);
               Hn(:, m) = Hn(:, m) + g * exp(-1j*2*pi*f*r0(l)/c) * an.';
               if m <= M
                   Gn(l, m) = g * exp(-1j*2*pi*f*r0(l)/c);
               end
=======
           if l == Ln
                g =  c / f / 4 / pi / r0(l);
                Hn(:, m) = Hn(:, m) + g * exp(-1j*2*pi*f*r0(l)/c) * an.';
                if m <= M
                   Gn(l, m) = g * exp(-1j*2*pi*f*r0(l)/c);
                end
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3
           %NLOS
           else
               g =  gn(l) * (c / f)^2 / (4 * pi)^2 / r0(l) / rsn(l);
               Hn(:, m) = Hn(:, m) + g * exp(-1j*2*pi*f*(r0(l)+rsn(l))/c) * an.';
               if m <= M
<<<<<<< HEAD
                   Gn(l, m) = g * exp(-1j*2*pi*f*(r0(l)+rsn(l))/c);
               end
           end 
        end
=======
                    Gn(l, m) = g * exp(-1j*2*pi*f*(r0(l)+rsn(l))/c);
               end
           end
       end 
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3
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
<<<<<<< HEAD

=======
%if Lf > 0
%    hf = hf * sqrt(Nt/Lf);
%end
%if Ln > 0
%    hn = hn * sqrt(Nt/Ln);
%end
%Cf=isequal(hf,h);
%Cn=isequal(hn,h);
>>>>>>> 316069b84f3774ffae0276fa8d15ee8c78e6faa3
end
