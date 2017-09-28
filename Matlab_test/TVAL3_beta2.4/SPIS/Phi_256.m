function [ W ] = Phi_256(ratio)
%Phi_256 return sensing matrix for 256 images
% Returns w/ <90% rows, Row_perm done w/ correct perm w/ correct col_perm

N = 2^16;
n = ((1:N)-1).';
gray = bin2gray(n,'pam',N);
clearvars n
br = de2bi(gray,'left-msb');
clearvars gray
s = bi2de(br) + 1;
s = single(s);
clearvars br
%

H = int8(hadamard(N/2));
%
H2 = zeros(N,N,'int8');
H2 = [H, H; H -H];
clearvars H
%

W = zeros(N,N,'int8');
W = H2(s,:);
clearvars s H2

%
W = uint8(W + 1);
%
W = W/2;
%
%Remove 40% due to max limit imwrite size
load row_perm_256.mat
load col_perm_256.mat
col_perm_wh = zeros(1,N);
col_perm_wh(col_perm) = 1:N;

if ratio > 0.9 || ratio < 0.05
    ratio = 0.9
end
M = round(N*ratio/24)*24;
picks = row_perm(1:M);
W = W(picks,:);
W = W(:,col_perm_wh);


end

