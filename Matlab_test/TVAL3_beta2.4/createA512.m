clear all
N = 2^18

%  'qam', 'pam', 'fsk', 'dpsk' or 'psk'
n = ((1:N)-1).';
gray = bin2gray(n,'pam',N);
clearvars n
br = de2bi(gray,'left-msb');
clearvars gray
s = bi2de(br) + 1;
s = single(s);
clearvars br
%

H128 = int8(hadamard(256^2/2));
%
H256 = zeros(256^2,256^2,'int8');
%
H256 = [H128, H128; H128 -H128];
%clearvars H128
%

H512 = zeros(N/16-2^6,N,'int8');
%
load 'row_perm_512.mat'

%
%perm s before constructing final matrix
picks = row_perm(1:(N/16-2^6));
s2 = s(picks);
%
for i = 1:(N/16-2^6)
    row = s2(i);
    in = mod(row,256^2) + 1;
    if row <= N/4
        H512(i,:) = [H256(in,:) H256(in,:) H256(in,:) H256(in,:)];
    elseif row <= N/2
        H512(i,:) = [H256(in,:) -H256(in,:) H256(in,:) -H256(in,:)];
    elseif row <= ((N/4)*3)
        H512(i,:) = [H256(in,:) H256(in,:) -H256(in,:) -H256(in,:)];
    else
        H512(i,:) = [H256(in,:) -H256(in,:) -H256(in,:) H256(in,:)];
    end
end


%
H512 = uint8(H512 + 1);
H512 = H512/2;
%%

imwrite(H512,'Psi/512_bin_row_perm_WH_10_2.PNG');