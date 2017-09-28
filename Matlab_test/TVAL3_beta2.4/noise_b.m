function [ bn ] = noise_b( b, n )
%return matrix w/ different amount of noise

N = length(b);
bn = zeros(N,5);
bn(:,1) = b;
% add noise
sigma = 0.00001;  % noise std

bavg = mean(abs(b));
noise = randn(N,1);

for m = 1:n
    bn(:,m+1) = b + sigma*bavg*noise*10^(m);
end

end

