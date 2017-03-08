% kc_example.m
%
% tests kc.m
%
% Created: Aug 2009.
% email: chinmay@rice.edu

clear, close all, clc
path(path,'../Utils')

N= 100;
K= 10;
C= 2;
iter = 10;
M = 26;

% generate signal
x = zeros(N,1); 
ind1= 10+(2:6);
ind2= 60+(10:14);
x(ind1,1)= (10+rand(length(ind1),1)*5).*sign(rand(length(ind1),1)-.5);
x(ind2,1)= (8+rand(length(ind2),1)*2).*sign(rand(length(ind2),1)-.5);

% measurements
Phi = (1/sqrt(M))*randn(M,N);
y = Phi*x(:);

% reconstruct
[xhat, trsh] = kc(y, Phi, K, C, iter);

figure(1), hold on
box on
stem(x,'.','MarkerSize',10)
stem(xhat,'rd','MarkerSize',3)
axisfortex('','','')
legend('Original','Reconstruction')