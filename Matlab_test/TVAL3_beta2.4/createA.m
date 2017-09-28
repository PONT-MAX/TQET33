clear all
N = 2^16
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
ratio = 0.6;
M = round(N*ratio);
picks = row_perm(1:M);
W = W(picks,:);
%%

imwrite(W,'Psi/256_bin_row_perm_WH_60_2.PNG');