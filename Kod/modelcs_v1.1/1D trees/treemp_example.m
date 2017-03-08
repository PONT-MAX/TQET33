% treemp_example.m
%
% tests treemp.m and treemp_fun.m
%
% IMPORTANT: Make sure you have the Rice Wavelet Toolbox in your path
%
% Created: Aug 2009.
% email: chinmay@rice.edu
clear, clc, close all
% Make signal
path(path,'../Utils')
N = 2^10;
h = daubcqf(4,'min');
L = log2(N);
MaxOrder = 3;
bpoints = 5;
% generate piecewise polynomial signal
s = transpose(genPWpoly(N,MaxOrder,bpoints));
sigma = 0; % signal noise

% K = 30;
iter = 20;
% wavelet decomposition

sn = s + randn(size(s))*sigma;
[wn0,lv] = mdwt(sn,h,L);
% figure(2), hold on, plot(s,'r', 'LineWidth', 3)

% tree-ify signal before feeding in
Bn = abs(wn0).^2;
maskn = cssa1(Bn,L,0.04*N);
wn=0*wn0; f=find(maskn>0);  wn(f) = wn0(f);
% wn = midwt(wC,h,L);

K = length(f); % size of tree
% M = round(7*K/2)*2;
M = 2*round(3/2*K);
% figure, plot(midwt(wC,h,L),'r')


% noisy measurements
Phi = (1/sqrt(M))*randn(M,N);
% y = PhiPsi(wn);
sg = 0.0;
%y = Phi*wn(:);
y = Phi*wn(:);
y = y + randn(size(y))*sg;

%----Reconstruction----%

% transform parameters
wvopts.filt = h; wvopts.lv = L; wvopts.length = N;
[xhat,trsh] = treemp(y, Phi, wvopts, K, iter);

figure(1), hold on
plot(midwt(wn,h,L),'r','LineWidth',2)
plot(xhat, 'b','LineWidth',3)

% Same as above, except that the operator \Phi is
% implemented using a function handle
% Use this for large signal sizes (N > 10000)
P = (1:N)';
q = randperm(N/2-1)+1; % this keeps a random set of FFT coefficients
OMEGA = q(1:M/2)';
Phi_f = @(z) A_f(z, OMEGA, P);
PhiT_f = @(z) At_f(z, N, OMEGA, P);

y2 = Phi_f(wn(:));

[xhat2,trsh] = treemp_fun(y2, Phi_f, PhiT_f, wvopts, K, iter);

figure(2), hold on
plot(midwt(wn,h,L),'r','LineWidth',2)
plot(xhat2,'b','LineWidth',3)