% jsmp_example.m
%
% tests jsmp.m and jsmp_fun.m
%
% Created: Aug 2009.
% email: chinmay@rice.edu

clear, close all, clc
path(path,'../Utils')

N = 4096; % signal length
J = 16; % signal sparsity
num_blocks = round(N/J);
K = round(0.05*num_blocks);
% so that overall sparsity = JK ~ 208

M = 600; % number of measurements

iter = 20;

% generate signal
v = randperm(num_blocks);
hi_col = v(1:K);
smat = zeros(J, num_blocks);
smat_hi = smat(:,hi_col);
smat_hi = randn(size(smat_hi));
smat(:,hi_col) = smat_hi;
x = smat(:);

% measurements
Phi = (1/sqrt(M))*randn(M,N);
y = Phi*x(:);

% reconstruct
[xhat, trsh] = jsmp(y,Phi,J,K,iter);

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

y2 = Phi_f(x(:));

[xhat2,trsh2] = jsmp_fun(y2, Phi_f, PhiT_f, N, K, J, iter);
figure(2), hold on
box on
stem(x,'.','MarkerSize',10)
stem(xhat2,'rd','MarkerSize',3)
axisfortex('','','')
legend('Original','Reconstruction')