% delta_example.m
%
% tests deltarec.m
%
% Created: Aug 2009.
% email:   chinmay@rice.edu

clear, close all, clc
path(path,'../Utils')

N= 1024;
K = 30; %sparsity
W = 30; % window size
iter = 30;

% generate test (K,W)-signal
% basically, generate a random sparse signal
% and perform a (K,W) approximation

w = randn(N,1); s = 0*w; 
v = randperm(N); s(v(1:K))=1; w=w.*s;
supp = bestdelta(w,K,W,1);
x = w.*round(supp);
% figure(123), stem(x_hat), pause
Ktru = sum(abs(x)>0)
M = round(4.5*Ktru);

% measurements
Phi = (1/sqrt(M))*randn(M,N);
y = Phi*x(:);

% reconstruct
[xhat, trsh] = deltarec(y, Phi, W, Ktru, iter);

figure(1), hold on
box on
stem(x,'.','MarkerSize',10)
stem(xhat,'rd','MarkerSize',3)
axisfortex('','','')
legend('Original','Reconstruction')