%% Init
addpath(genpath('/Users/Andreas/Skola/TQET33/Kod'))
%% 
% Vipe
clear, clc, close all, clear all

names = {'disguise.jpg', 'driver.jpg', 'fire.jpg', 'glass.jpg',...
    'industrial.jpg', 'lowlight.jpg', 'maritime.jpg',...
    'produce.jpg', 'rain.jpg', 'tornadoresearch.jpg'};

path(path,'../Utils')
filename = names{1};
s = imresize(double(rgb2gray(imread(filename))),0.25, 'bicubic');
s = s./max(max(s));
%figure, imagesc(s); colormap(gray(256)); axis image, rmaxis, shg

[n,m] = size(s);
N = m*n;
M = round(N/4);
K = round(M/5);
iter = 50;

Phi = randi(2, [M,N]) - 1;

y = phi*s(:);
