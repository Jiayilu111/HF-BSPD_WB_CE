function  at = polar_domain_manifold( Nt, d, f, theta0, alpha0 )
    c = 3e8;
    nn = -(Nt-1)/2:1:(Nt-1)/2;
%     r = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sin(theta0));
%     r = r0 - nn * d * sin(theta0) + nn.^2 * d^2 * cos(theta0)^2 ./ 2 /r0;
    r =  - nn * d * sin(theta0) + nn.^2 * d^2 * alpha0;
    at = exp(-1j*2*pi*f*(r)/c)/sqrt(Nt);
end

