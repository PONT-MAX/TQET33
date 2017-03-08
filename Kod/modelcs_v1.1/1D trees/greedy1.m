
% GREEDY1.m:  GREEDY WAVELET TREE ALGORITHM (1-d)
%
% B:     input data
% L:     number of levels in wavelet transform
% vol:   stopping volume of mask (default is full volume)
%        NOTE: does NOT include scaling coeff volume!
%
% mask:  final greedy energy shell configurations, ordered 
%
% Assumes that the underlying tree is the 1-d dyadic wavelet tree
%
% RGB INI September 1998

function mask = greedy1(B,L,vol);

%---------------------------------------------------------------------------%
% INITIALIZATIONS

% make input data a column vector
B    = B(:);
N    = length(B);
cvol = 0;
sno  = 2;
% scaling coefficient with highest index in wavelet array
parroot = N/2^L;

if nargin < 3
  vol = (N-N/2^L);
end
if vol > (N-N/2^L)
  disp('ERROR: volume parameter cannot be > number of wavelet coeffs')
  mask = [];
  return
end


% BQ: array with B vals and wc positions
BQ(1:N,1) = (1:N)';
BQ(1:N,2) = B;
BQ(1:N/2^L,2) = NaN*ones(N/2^L,1);

% Set all scaling coeffs to have mask = 1
mask = zeros(N,1);               % output mask is zero at start
mask(1:N/2^L) = ones(N/2^L,1)';  % except for parents of roots ("zero" nodes)

% don't exactly know why we need this, since we should always have vv=0
vv = length(find(mask(N/2^L+1:N)>0));

%---------------------------------------------------------------------------%
% MAIN LOOP (follow pseudo code on p. 138 and steps on p. 142)

while and( length(find(mask==0))>0 , vv<vol )

  % STEP 1
  % find node with largest value
  % delete it from the BQ as well
  maxn = Qfind(BQ);
  BQ   = Qdel(BQ,maxn);
  
  % STEP 2
  % form greedy supernode by traveling from parent to parent until
  % we hit a gsn where mask is not zero
  mask(maxn) = sno;
  par = ceil(maxn/2);
  while mask(par) == 0
    mask(par) = sno;
    BQ = Qdel(BQ,par);
    par = ceil(par/2);
  end
  sno = sno + 1;
  
  vv = length(find(mask(N/2^L+1:N)>0));
end
return