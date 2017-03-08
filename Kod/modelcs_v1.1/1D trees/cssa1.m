
% CSSA1.m:  CONDENSING SORT AND SELECT ALGORITHM (1-d)
%
% B:     input data
% L:     number of levels in wavelet transform
% vol:   stopping volume of mask (default is full volume)
%        NOTE: does NOT include scaling coeff volume!
%
% mask:  final supernode configurations, ordered 
%
% Assumes that the underlying tree is the 1-d dyadic wavelet tree
%
% RGB INI September 1998

function mask = cssa1(B,L,vol);

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

% set up each node as a supernode.  use an Nx5 matrix:
%    col 1: utp  (uptree pointer)
%    col 2: rp   (supernode root pointer)
%    col 3: snv  (supernode value)
%    col 4: num  (number of internal nodes)
%    col 5: mask (mask/kernel value)

sn(1:N, 1) = (1:N)';      % point uptree top pointers to self
sn(1:N, 2) = -ones(N,1);  % each node is a sn root => -1 flag
sn(1:N, 3) = B;           % each snv is the data itself
sn(1:N, 4) = ones(N,1);   % num = 1 for each sn   
sn(1:N, 5) = zeros(N,1);  % output mask is zero at start
sn(1:N/2^L, 5) = ones(N/2^L,1);  % except for parents of roots ("zero" nodes)

% set up SNV queue  (include scaling coeff values just so that Q never empties
snvQ(1:N, 1) = (1:N)';    % each node is its own sn topnode
snvQ(1:N, 2) = B;         % each snv is the wavelet data itself
snvQ(1:N/2^L,2) = NaN*snvQ(1:N/2^L,2);


%---------------------------------------------------------------------------%
% MAIN LOOP (follow pseudo code on p. 138 and steps on p. 142)

% notation:   MAX refers to the current sn we take from snv queue
%             CHECK refers to the parent sn it abutts

while (cvol < vol)

  % STEP 1
  % find supernode with largest SNV (sn MAX)
  % delete it from the snvQ as well
  topmax = Qfind(snvQ);
  snvQ = Qdel(snvQ,topmax);
  
  % STEP 2
  % find the root node of sn MAX and then its parent
  % this latter node is a leaf of the CHECK sn above ours
  if sn(topmax,2) == -1
    % then topmax is the sn root
    leafcheck = ceil(topmax/2);
  else
    % go to the sn root
    leafcheck = ceil(sn(topmax,2)/2);
  end
  
  
  % check whether CHECK sn has mask=1,0 (first find top node of CHECK)
%  topcheck = Tfind(sn,leafcheck);
  [topcheck,sn] = Tfind(sn,leafcheck);       % need to pass sn back
                                             % if we do path comp step
    
  if sn(topcheck,5) > 0 
    % we are in the money - CHECK sn has mask > 0.  so take sn MAX
    % STEP 3 (usual case - not the first sn taken): 
    sn = umask(sn,topmax,sno,N);
    cvol = cvol + sn(topmax,4); 
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

mask = sn(:,5);
return