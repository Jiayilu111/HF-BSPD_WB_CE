function hn = channel_norm(h)
[N, M] = size(h);

% hn = zeros(size(h));

hn = h / sqrt(sum(abs(h(:)).^2)) * sqrt(N * M);


end

