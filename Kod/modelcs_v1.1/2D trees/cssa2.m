
% CSSA2.m:  CONDENSING SORT AND SELECT ALGORITHM (2-d)
%
% B:     input data matrix
% L:     number of levels in wavelet transform
% vol:   stopping volume of mask (default is full volume)
%        NOTE: does NOT include scaling coeff volume!
%
% mask:  final supernode configurations, ordered 
%
% Assumes that the underlying tree is the 2-d dyadic wavelet tree
%
% RGB Rice January 1999

function mask = cssa2(B,L,vol);

%---------------------------------------------------------------------------%
% INITIALIZATIONS

[M,N] = size(B);  % assume that M=N (square matrix B)
cvol = 0;
sno  = 2;
% row and column index of scaling coefficient with highest index 
%  in wavelet array
parroot = N/2^L;

if nargin < 3
  vol = N^2 - (N/2^L)^2 ;
end
if vol > (N^2 - (N/2^L)^2)
  disp('ERROR: volume parameter cannot be > number of wavelet coeffs')
  mask = [];
  return
end

% set up each node as a supernode.  use an NxNx7 matrix:
%    last col 1: utpR  (uptree row pointer)
%    last col 2: utpC  (uptree column pointer)
%    last col 3: rpR   (supernode root row pointer)
%    last col 4: rpC   (supernode root column pointer)
%    last col 5: snv   (supernode value)
%    last col 6: num   (number of internal nodes)
%    last col 7: mask  (mask/kernel value)

sn = zeros(N,N,7);
sn(:,:,1) = (1:N)' * ones(1,N);  % point uptree top R pointers to self
sn(:,:,2) = ones(N,1) * (1:N);   % point uptree top C pointers to self
sn(:,:,3) = -ones(N,N);  % each node is a sn root => -1 flag
sn(:,:,4) = -ones(N,N);  % each node is a sn root => -1 flag
sn(:,:,5) = B;             % each snv is the data itself
sn(:,:,6) = ones(N,N);   % num = 1 for each sn   
sn(:,:,7) = zeros(N,N);  % output mask is zero at start
sn(1:N/2^L, 1:N/2^L, 7) = ones(N/2^L, N/2^L);   % except for 
                                 % parents of roots ("zero" nodes)

% set up SNV queue as an N^2 x 3 array
% (include scaling coeff values just so that Q never empties; make
%  then NaNs)
%    col 1: utpR  (uptree row pointer)
%    col 2: utpC  (uptree column pointer)
%    col 3: snv   (supernode value)
Bt = B; 
Bt(1:N/2^L,1:N/2^L) = NaN*ones(N/2^L, N/2^L);
r = sn(1:N, 1:N, 1);  c = sn(1:N, 1:N, 2);   % each node is its own sn topnode
snvQ = [r(:), c(:), Bt(:)];


%---------------------------------------------------------------------------%
% MAIN LOOP (follow pseudo code on p. 138 and steps on p. 142)

% notation:   MAX refers to the current sn we take from snv queue
%             CHECK refers to the parent sn it abutts
kk = 1
while (cvol < vol)
%   kk = kk+1, size(snvQ)
  % STEP 1
  % find supernode with largest SNV (sn MAX)
  % delete it from the snvQ as well
  %size(snvQ)
  topmax = Qfind(snvQ);
  topmaxR = topmax(1); topmaxC = topmax(2);
  snvQ = Qdel(snvQ,topmax);
%   size(snvQ)
  % STEP 2
  % find the root node of sn MAX and then its parent
  % this latter node is a leaf of the CHECK sn above ours
  if sn(topmaxR,topmaxC,3) == -1
    % then topmax is the sn root
    leafcheckR = ceil(topmaxR/2);
    leafcheckC = ceil(topmaxC/2);
  else
    % go to the sn root
    leafcheckR = ceil(sn(topmaxR,topmaxC,3)/2);
    leafcheckC = ceil(sn(topmaxR,topmaxC,4)/2);
  end
  leafcheck = [leafcheckR, leafcheckC];

  % check whether CHECK sn has mask=1,0 (first find top node of CHECK)
  [topcheck, sn] = Tfind(sn,leafcheck);
  topcheckR = topcheck(1); topcheckC = topcheck(2); 
  
  if sn(topcheckR,topcheckC,7) > 0 
    % we are in the money - CHECK sn has mask > 0.  so take sn MAX
    % STEP 3 (usual case - not the first sn taken): 
    sn = umask(sn,topmax,sno,N);
        
    cvol = cvol + sn(topmaxR,topmaxC,6); 
    sno = sno + 1;
 
  else
    % CHECK sn has mask = 0, so merge CHECK and MAX sn's
    % STEP 4
    [sn,snvQ] = condense(sn,snvQ,topcheck,topmax);
    
  end
%  sn
%  snvQ
%  pause
end

mask = sn(:,:,7);
return
