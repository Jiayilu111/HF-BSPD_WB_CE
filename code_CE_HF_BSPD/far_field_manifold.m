function  at = far_field_manifold(Nt,theta,d, lambda)
    at = exp(1i*2*pi*[0:Nt-1]*d*sin(theta)/lambda);
    at = at / sqrt(Nt);
end
