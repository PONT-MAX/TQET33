
% GREEDY2.m:  GREEDY WAVELET TREE ALGORITHM (2-d)
%
% B:     input data
% L:     number of levels in wavelet transform
% vol:   stopping volume of mask (default is full volume)
%        NOTE: does NOT include scaling coeff volume!
%
% mask:  final greedy energy shell configurations, ordered 
%
% Assumes that the underlying tree is the 2-d dyadic wavelet tree
%
% RGB INI September 1998
% RGB Rice  January 1999

function mask = greedy2(B,L,vol);

%---------------------------------------------------------------------------%
% INITIALIZATIONS

[M,N] = size(B);  % assume that M=N (square matrix B)
cvol = 0;
sno  = 2;
% scaling coefficient with highest index in wavelet array
parroot = N/2^L;

if nargin < 3
  vol = N^2 - (N/2^L)^2 ;
end
if vol > (N^2 - (N/2^L)^2)
  disp('ERROR: volume parameter cannot be > number of wavelet coeffs')
  mask = [];
  return
end

% set up BQ queue as an N^2 x 3 array
% (include scaling coeff values just so that Q never empties; make
%  then NaNs)
%    col 1: R  (row pointer)
%    col 2: C  (column pointer)
%    col 3: B  (B value)
Bt = B; 
Bt(1:N/2^L,1:N/2^L) = NaN*ones(N/2^L, N/2^L);
r = (1:N)' * ones(1,N);  
c = ones(N,1) * (1:N);   			
BQ = [r(:), c(:), Bt(:)];

% Set all scaling coeffs to have mask = 1
mask = zeros(N,N);               % output mask is zero at start
mask(1:N/2^L,1:N/2^L) = ones(N/2^L,N/2^L);  % except for parents 
                                 % of roots ("zero" nodes)

% don't exactly know why we need this, since we should always have vv=0
%vv = length(find(mask(N/2^L+1:N,N/2^L+1:N)>0));
vv = 0;

%---------------------------------------------------------------------------%
% MAIN LOOP

while and( length(find(mask==0))>0 , vv<vol )

  % STEP 1
  % find node with largest value
  % delete it from the BQ as well
  maxn = Qfind(BQ);
  maxnR = maxn(1); maxnC = maxn(2);
  BQ   = Qdel(BQ,maxn);
  
  % STEP 2
  % form greedy supernode by traveling from parent to parent until
  % we hit a gsn where mask is not zero
  mask(maxnR,maxnC) = sno;
  parR = ceil(maxnR/2);
  parC = ceil(maxnC/2);
  par  = [parR, parC];
 
  while mask(parR,parC) == 0
    mask(parR,parC) = sno;
    BQ = Qdel(BQ,par);
    parR = ceil(parR/2);
    parC = ceil(parC/2);
    par = [parR, parC];
  end
  sno = sno + 1;
  
  % volume of mask (excluding scaling coeffs)
  vv = length(find(mask>0)) - (N/2^L)^2;
end
return