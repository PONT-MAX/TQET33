
% Qswap.m : Delete one element and insert another into the SNV queue 
%           (1-d or 2-d)
%
% snvQout = Qswap(snvQin,topnodeD,snvI,topnodeI)
%
% snvQin:    SNV queue 
%            1-d: top node in column 1, snv in column 2
%            2-d: top node row in column 1, column in column 2, snv in column 3
% topnodeD:  topnode of the element to be deleted
% snvI:      SNV of the element to be inserted
% topnodeI:  topnode of the element to be inserted
%
% snvQout:  updated SNV queue
%
% This command is useful only because in the current implementation
% following Qdel.m by Qins.m is inefficient
%
% RGB INI September 1998

function snvQout = Qswap(snvQin,topnodeD,topnodeI,snvI);

[N,d] = size(snvQin);
topnodeD = topnodeD(:);
topnodeI = topnodeI(:);

if d==2
  snvQout = snvQin;
  % find topnodeD
  here  = find(snvQin(:,1)==topnodeD);
  if isempty(here)
    disp('ERROR: Qswap did not find a supernode with this topnodeD')
    disp('topnodeD'); topnodeD
    snvQout = [];
    return
  end
  % swap
  snvQout(here,:) = [topnodeI,snvI];
  
else
  
  snvQout = snvQin;
  % find topnodeD (thanks to jrom!)
  hereM = (snvQin(:,1:2)==repmat(topnodeD',[N 1]));
  here = find((hereM(:,1).*hereM(:,2))==1);
   
  if isempty(here)     % not tested in 2-d
    disp('ERROR: Qswap did not find a supernode with this topnode')
    disp('topnodeD'); topnodeD
    snvQout = [];
    return
  end
  
  % swap
  snvQout(here,:) = [topnodeI',snvI];
  
end
return



