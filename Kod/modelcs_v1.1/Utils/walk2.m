
% WALK2.m:  Step one scale through a 2-d wavelet sn subtree starting
%             at the root and set mask = 1 at each node
%           When called recursively, will look thru the entire
%             wavelet sn subtree
%
% snout = walk2(snin,root,N);
%
% snin:   supernode data base
% root:   current root node (explore from here down) (row and column)
% sno:    number to put on mask (orders supernodes)
% N:      max index in wavelet array
%
% snout:  updated supernode data base
%
% ASSUMES wavelet coeffecients are numbered 1:N on rows/columns
%
% RGB Rice January 1999

function snout = walk2(snin,root,sno,N);

sn = snin;
rootR = root(1); rootC = root(2);

lchildR = 2*rootR - 1;
rchildR = 2*rootR;

lchildC = 2*rootC - 1;
rchildC = 2*rootC;

% THIS STEP COULD BE WRONG IN 2-D!!!!!!!! CHECK IT!
if or(lchildR>N, lchildC>N)  % stop traversal when we get to leaves of
  snout = sn;                % wavelet tree
%  disp('bailing out'), lchildR, lchildC
  return
end

if sn(lchildR,lchildC,3) ~= -1        % -1 in sn root field means sn root
  rootR = lchildR;  rootC = lchildC;  root = [rootR, rootC];
  sn(rootR,rootC,7) = sno;
  sn = walk2(sn,root,sno,N);     % keep walking down
end

if sn(rchildR,lchildC,3) ~= -1        % -1 in sn root field means sn root
  rootR = rchildR;  rootC = lchildC;  root = [rootR, rootC];
  sn(rootR,rootC,7) = sno;
  sn = walk2(sn,root,sno,N);     % keep walking down
end

if sn(lchildR,rchildC,3) ~= -1        % -1 in sn root field means sn root
  rootR = lchildR;  rootC = rchildC;  root = [rootR, rootC];
  sn(rootR,rootC,7) = sno;
  sn = walk2(sn,root,sno,N);     % keep walking down
end

if sn(rchildR,rchildC,3) ~= -1        % -1 in sn root field means sn root
  rootR = rchildR;  rootC = rchildC;  root = [rootR, rootC];
  sn(rootR,rootC,7) = sno;
  sn = walk2(sn,root,sno,N);     % keep walking down
end

snout = sn;
return



