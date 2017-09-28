%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A_fWH
%
% Written by: Chengbo Li
% CAAM, Rice University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = A_fWH(x, OMEGA, P)

% A_fWH(x,picks,perm);

N = length(x);

x = x(:);
fx = fWHtrans(x(P))*sqrt(N);
b = fx(OMEGA);