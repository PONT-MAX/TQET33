
% UMASK.m:  Update mask field in a taken supernode (to all 1's)  (1-d or 2-d)
%
% snout = umask(snin,topnode,sno,N);
%
% snin:      supernode data base
% topnode:   top node of supernode
% sno:       number to put on mask (orders supernodes)
% N:         max index in wavelet array
%
% snout:     updated supernode data base
%
% RGB INI September 1998

function snout = umask(snin,topnode,sno,N);

if length(topnode)==1    % 1-d case
  sn = snin;
  % traverse sn subtree starting from root using recursive function
  % set mask=sno at each node
  if sn(topnode,2) ~= -1           % root and top are different nodes
    root = sn(topnode,2);
  else                             % root and top are the same node
    root = topnode;
  end
  sn(root,5) = sno;
  sn         = walk1(sn,root,sno,N);
  snout = sn;
  
else                               % 2-d case
  sn = snin;
  topR = topnode(1); topC = topnode(2);
  % traverse sn subtree starting from root using recursive function
  % set mask=sno at each node
  if sn(topR,topC,3) ~= -1         % root and top are different nodes
    rootR = sn(topR,topC,3);
    rootC = sn(topR,topC,4);
  else                             % root and top are the same node
    rootR = topR;  rootC = topC;
  end
  sn(rootR,rootC,7) = sno;
  root = [rootR, rootC];
  sn = walk2(sn,root,sno,N);
  snout = sn;
 
end
return



