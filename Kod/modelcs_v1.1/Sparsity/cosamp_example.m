% cosamp_example.m
%
% tests cosamp.m and cosamp_fun.m
%
%
% Created: Aug 2009.
% email: chinmay@rice.edu

clear, close all, clc
path(path,'../Utils')

N = 1024; % signal length
K = 40; % signal sparsity
M = 240; % number of measurements

iter = 10; % number of Cosamp iterations

% generate signal
x = zeros(N,1);
v = randperm(N);
x(v(1:K)) = randn(K,1);

% measurements
Phi = (1/sqrt(M))*randn(M,N);
y = Phi*x;

% reconstruct
[xhat, trsh] = cosamp(y,Phi,K,iter);

figure(1), hold on
box on
stem(x,'.','MarkerSize',10)
stem(xhat,'rd','MarkerSize',3)
axisfortex('','','')
legend('Original','Reconstruction')

% Same as above, except that the operator \Phi is
% implemented using a function handle
% Use this for large signal sizes (N > 10000)
P = (1:N)';
q = randperm(N/2-1)+1; % this keeps a random set of FFT coefficients
OMEGA = q(1:M/2)';
Phi_f = @(z) A_f(z, OMEGA, P);
PhiT_f = @(z) At_f(z, N, OMEGA, P);

y2 = Phi_f(x);

[xhat2,trsh] = cosamp_fun(y2, Phi_f, PhiT_f, N, K, iter);
figure(2), hold on
box on
stem(x,'.','MarkerSize',10)
stem(xhat2,'rd','MarkerSize',3)
axisfortex('','','')
legend('Original','Reconstruction')