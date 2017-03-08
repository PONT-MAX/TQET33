% script to test tree based recovery on images
% 2D version
% Chin, Aug 13 '08

% IMPORTANT: make sure you have the Rice Wavelet Toolbox in your path
% http://dsp.rice.edu/software/rice-wavelet-toolbox

% reistributed: August '09
% slow and buggy. Use with care :(
% email: chinmay@rice.edu
%%
addpath(genpath('/Users/Andreas/Skola/TQET33/Kod'))
%%
clear, clc, close all, clear all
% load image
path(path,'../Utils')
names = {'disguise.jpg', 'driver.jpg', 'fire.jpg', 'glass.jpg',...
    'industrial.jpg', 'lowlight.jpg', 'maritime2.jpg',...
    'produce.jpg', 'rain.jpg', 'tornadoresearch.jpg'};
scol = imread(['../../Ex_bilder/',names{7}]);
sg = imresize(double(rgb2gray(scol)),0.25, 'bicubic');
s = sg./max(max(sg));

% compute wavelet transform
[n,m] = size(s);
sigma = 0;
N = m*n;
% K = 30;
iter = 50;
% wavelet decomposition
L=6;
h = daubcqf(8,'min');
sn = s + randn(size(s))*sigma;
[wn0,lv] = mdwt(sn,h,L);
% figure(2), hold on, plot(s,'r', 'LineWidth', 3)

% tree-ify signal before feeding in
Bn = abs(wn0).^2;
maskn = greedy2(Bn,L,.10*N);
wn=0*wn0; f=find(maskn>0);  wn(f) = wn0(f);
[snew, lv] = midwt(wn, h, L);

K = length(f); % size of tree
% M = round(7*K/2)*2;
M = round(4*K);
% figure, plot(midwt(wC,h,L),'r')

% noisy measurements
%Phi = randi(2,[M,N]) - 1;
Phi = (1/sqrt(M))*randn(M,N);
% y = PhiPsi(wn);
sg = 0*0.1;
%y = Phi*wn(:);
y = Phi*wn0(:);
y = y + randn(size(y))*sg;

%Bn = abs(wn).^2;
%maskn = cssa1(Bn,L,vol);

%----Reconstruction----%
%
% transform parameters
wvopts.filt = h; wvopts.lv = L; wvopts.xlength = n;
wvopts.slength = N; % assume square 


%%%%%%%%%%%%%%%%%%%%%
% Note: 
% treemp_nf_2D.m and treemp_greedy_2D.m both do tree-approximations of 
% sparse signals
% treemp_nf_2D uses CSSA which is theoretically optimal, but is buggy
% If it fails, use treemp_greedy_2D
%%%%%%%%%%%%%%%%%%%%%


xhat_nf = treemp_nf_2D(y, Phi, wvopts, K, iter, wn);

xhat = treemp_greedy_2D(y, Phi, wvopts, K, iter, wn0(:));

% /Note

norm(xhat(:,end)-sn(:))
xhat2 = cosamp_2D_nf(y, Phi, wvopts, K, iter, wn0(:));
norm(xhat2(:,end)-sn(:))
%figure, imagesc(snew); colormap(gray(256)); axis image, rmaxis, shg
figure(4), imagesc(s); colormap(gray(256)); axis image, rmaxis, shg