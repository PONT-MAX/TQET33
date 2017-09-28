%% How to get the same transfor that TVAL3 uses in its fast mode.

side_length = 8;
N = side_length^2;
%load('row_perm_128');
%load('col_perm_128');
% generate measurement matrix
p = randperm(N); % Generera stor matris
ratio = 1;
M = round(N*ratio);
picks = p(1:M); % Psi välj M av p
% En kontroll för olämplia val?
%for ii = 1:M
%    if picks(ii) == 1 % cant have row 1 in hw-matrix
%        picks(ii) = p(M+1);
%        break;
%    end
%end
clearvars p
%
col_perm = randperm(N); % column permutations allowable
col_perm_wh = zeros(1,N);
col_perm_wh(col_perm) = 1:N;

A = @(x,mode) dfA(x,picks,col_perm,mode);
Im = eye(side_length);

% Create WH-basis matrix
WH = fWHtrans(eye(N)*sqrt(N));

%WH2 = A(I,1);
%%
Im2 = Im(:);
% W/ matri1x 11mult
B = WH(:,col_perm_wh)*Im2;
B = B(picks);
% TVAL3 method
b = A(Im2,1);

% Show res
[b(1:6) B(1:6)]
isequal(b,B)

%%
clc
perm = [2 3 1 4 5]
a = [11 33 55 17 19]'
a(perm)
%%
perm2 = res;
b = [1 1 0 1 1; 0 0 1 0 0; 0 1 0 1 0]
b(:,perm2)
%%
b*a
b*a(perm)
b(:,perm2)*a

%%

N = length(perm);
n = 1:N;
s = perm - n;
res = zeros(1,N);
res2 = zeros(1,N);
res2(perm) = 1:N
for i = 1:N
    res(i + s(i)) = i;
end
perm
res


