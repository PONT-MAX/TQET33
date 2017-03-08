
% WALK1.m:  Step one scale through a 1-d wavelet sn subtree starting
%             at the root and set mask = 1 at each node
%           When called recursively, will look thru the entire
%             wavelet sn subtree
%
% snout = walk1(snin,root,N);
%
% snin:   supernode data base
% root:   current root node (explore from here down)
% sno:    number to put on mask (orders supernodes)
% N:      max index in wavelet array
%
% snout:  updated supernode data base
%
% ASSUMES wavelet coeffecients are numbered 1:N
%
% RGB INI September 1998

function snout = walk1(snin,root,sno,N);

sn = snin;
lchild = 2*root - 1;
rchild = 2*root;

if lchild > N                % stop travesal when we get to leaves of
  snout = sn;                % wavelet tree
%  disp('bailing out'), lchild
  return
end

if sn(lchild,2) ~= -1        % -1 in sn root field means sn root
  root = lchild;
  sn(root,5) = sno;
  sn = walk1(sn,root,sno,N);     % keep walking down
end

if sn(rchild,2) ~= -1        % -1 in sn root field means sn root
  root = rchild;
  sn(root,5) = sno;
  sn = walk1(sn,root,sno,N);     % keep walking down
end

snout = sn;
return



