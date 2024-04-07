Np = 1000; Nq = 1000;
p_list = linspace(-0.01, 0.01, Np);
q_list = linspace(-0.05, 0.05, Nq);

x = linspace(-0.25, 0.25, 1000);

tic()
y = zeros(Np, Nq);
for m = 1:Np
    p = p_list(m);
    for n = 1:Nq
        q = q_list(n);
        f = exp( 2j * pi / 0.001 * ( p * x - q * x.^2 ) );
        y(m, n) = sum(f) * (x(2) - x(1)) * 2; %Control the correlation function within 0-1
    end 
end
toc()

figure;
mesh(q_list, p_list, abs(y))
xlabel('$\eta_{l,c}$-$\frac{f_m}{f_c}\eta_{l,m}$', 'interpreter', 'latex')
ylabel('$\theta_{l,c}$-$\frac{f_m}{f_c}\theta_{l,m}$', 'interpreter', 'latex')
zlabel('$G_{near}$', 'interpreter', 'latex')
colormap('jet')
colorbar
