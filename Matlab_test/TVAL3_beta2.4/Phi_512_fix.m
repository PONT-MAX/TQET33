function [ H512 ] = Phi_512_fix( ratio, start, stop )
%Phi_512 return sensing matrix for 512 images
% Returns w/ <N^2/16% rows, Row_perm & Col_perm done w/ correct perm
N = 2^18;

% Calculates Index perm for Hadamard to Seq orderd WH-matrix
n = ((1:N)-1).';
% Gray code
gray = bin2gray(n,'pam',N);
clearvars n
% Reverse binary
br = de2bi(gray,'left-msb');
clearvars gray
s = bi2de(br) + 1;
s = single(s);
clearvars br
%

% Sart w/ small to not run out of mem
H128 = int8(hadamard(256^2/2));
%
H256 = zeros(256^2,256^2,'int8');
% Cunstruct 256^2
N1 = 256^2/2;
H256(1:N1,1:N1) = H128;
clearvars H128
%
H256(N1+1:end,1:N1) = H256(1:256^2/2,1:256^2/2);
H256(1:N1,N1+1:end) = H256(1:256^2/2,1:256^2/2);
H256(N1+1:end,N1+1:end) = -H256(1:256^2/2,1:256^2/2);
%
% This is the max limit before memory runs out
if ratio > 0.35 || ratio < 0
    ratio = 0.35
end
M = round(N*ratio/9)*9;
H512 = zeros(M/2,N,'int8');
%
load 'row_perm_512.mat'
load col_perm_512.mat
col_perm_wh = zeros(1,N);
col_perm_wh(col_perm) = 1:N;
clearvars col_perm
%
%perm s before constructing final matrix
picks = row_perm(1:M);
s2 = s(picks);
clearvars picks row_perm
ins = zeros(1,N,'int8');
%
k = 1;
for i = start:stop
    if mod(i,500) == 0
        100*double(i)/double(M)
    end
    row = s2(i);
    in = mod(row,256^2);
    if in == 0
        in = 256^2;
    end
    if row <= N/4
        ins = [H256(in,:) H256(in,:) H256(in,:) H256(in,:)];
        ins = ins(:,col_perm_wh);
    elseif row <= N/2
        ins = [H256(in,:) -H256(in,:) H256(in,:) -H256(in,:)];
        ins = ins(:,col_perm_wh);
    elseif row <= ((N/4)*3)
        ins = [H256(in,:) H256(in,:) -H256(in,:) -H256(in,:)];
        ins = ins(:,col_perm_wh);
    else
        ins = [H256(in,:) -H256(in,:) -H256(in,:) H256(in,:)];
        ins = ins(:,col_perm_wh);
    end
    
    H512(k,:) = ins;
    k = k + 1;
end
clearvars H256 
%
H512 = (H512 + 1);
H512 = H512/2;

end

