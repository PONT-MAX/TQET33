
% Qfind.m : Find the largest element on the SNV queue (1-d or 2-d)
%
% [topnode,snv] = Qfind(snvQ)
%
% snvQ:     SNV queue 
%           1-d: top node in column 1, snv in column 2
%           2-d: top node row in column 1, column in column 2, snv in column 3
%
% topnode:  top node of this supernode
% snv:      maximum SNV on the SNV queue
%
% RGB Rice January 1999

function [topnode,snv] = Qfind(snvQ);

[l,d] = size(snvQ);
if d==2
  [m,loc] = max(snvQ(:,2));
  topnode = snvQ(loc,1);
  snv     = snvQ(loc,2);

else
  [m,loc] = max(snvQ(:,3));
  topnode = [snvQ(loc,1),snvQ(loc,2)];
  snv     = snvQ(loc,3);

end
return



