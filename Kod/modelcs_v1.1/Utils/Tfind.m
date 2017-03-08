
% Tfind.m : Find the top node of a supernode (1-d or 2-d)
%
% topnode = Tfind(sn,node)
%
% sn:       supernode data base
% node:     node from supernode for which we wish to find top node
%           (a vector with row and column in 2-d case)
%
% topnode:  topnode of the supernode containing node
%           (a vector with row and column in 2-d case)
%
% RGB Rice January 1999

function [topnode,sn] = Tfind(sn,node);

if length(size(sn))==2         % 1-d case
  % simply jump from uptree pointer to uptree pointer 
  % until we point to ourselves (=> we are at the top node)
  up = sn(node,1);
  while sn(up,1) ~= up
    up = sn(up,1);
  end
  topnode = up;

  % path compression step (see p. 141):
  % now that we know where the topnode is, go thru path again setting
  % all up pointers to topnode
  up         = sn(node,1);
  sn(node,1) = topnode;
  while sn(up,1) ~= up
    up1 = sn(up,1);
    sn(up,1) = topnode;
    up = up1;
  end
   
elseif length(size(sn))==3      % 2-d case
  % simply jump from uptree pointer to uptree pointer 
  % until we point to ourselves (=> we are at the top node)
  nodeR = node(1); 
  nodeC = node(2);
  upR = sn(nodeR,nodeC,1); 
  upC = sn(nodeR,nodeC,2);
  
  while or(sn(upR,upC,1)~=upR, sn(upR,upC,2)~=upC)
    upRtmp = sn(upR,upC,1);
    upC    = sn(upR,upC,2);
    upR    = upRtmp;
  end
  topnode = [upR, upC];

  % path compression step (see p. 141):
  % now that we know where the topnode is, go thru path again setting
  % all up pointers to topnode
%  nodeR = node(1);                 
%  nodeC = node(2);
  upR = sn(nodeR,nodeC,1);          
  upC = sn(nodeR,nodeC,2);
  sn(nodeR,nodeC,1) = topnode(1);   
  sn(nodeR,nodeC,2) = topnode(2);
  
  while or(sn(upR,upC,1)~=upR, sn(upR,upC,2)~=upC)
    up1R = sn(upR,upC,1);           up1C = sn(upR,upC,2);
    sn(upR,upC,1) = topnode(1);     sn(upR,upC,2) = topnode(2);
    upR = up1R;                     upC = up1C;
  end
 
end
return



