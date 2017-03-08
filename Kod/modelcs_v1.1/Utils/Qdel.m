
% Qdel.m : Delete an element from the SNV queue (1-d or 2-d)
%
% snvQout = Qdel(snvQin,topnode)
%
% snvQin:   SNV queue
%           1-d: top node in column 1, snv in column 2
%           2-d: top node row in column 1, column in column 2, snv in column 3
%
% topnode:  topnode of the element to be deleted
%
% snvQout:  updated SNV queue
%
% RGB Rice January 1999

function snvQout = Qdel(snvQin,topnode);

[N,d] = size(snvQin);
if N<2
  disp('Qdel has deleted final entry from the queue')
  snvQout = []; return	
end
		      
if d==2             % 1-d
  % find topnode
  here  = find(snvQin(:,1)==topnode);

  if isempty(here)
    disp('ERROR: Qdel did not find a supernode with this topnode')
    snvQout = [];
    return
  end

  % move entries from here+1 to N up one place
  snvQout(1:here-1,:) = snvQin(1:here-1,:);
  if here < N
    snvQout(here:N-1,:) = snvQin(here+1:N,:);
  end
  
  
else               % 2-d
  % find topnode (thanks to jrom!)
  hereM = (snvQin(:,1:2)==repmat(topnode,[N 1]));
  here = find((hereM(:,1).*hereM(:,2))==1);

  if isempty(here)     % not tested in 2-d
    disp('ERROR: Qdel did not find a supernode with this topnode')
    topnode
    snvQout = [];
    return
  end

  % move entries from here+1 to N up one place
  snvQout(1:here-1,:) = snvQin(1:here-1,:);
  if here < N
    snvQout(here:N-1,:) = snvQin(here+1:N,:);
  end
  
end
return