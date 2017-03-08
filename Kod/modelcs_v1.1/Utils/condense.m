
% CONDENSE.m:  Condense two supernodes into one   (1-d or 2-d)
%
% [snout,snvQout] = condense(snin,snvQin,top1,top2);
%
% snin:      supernode data base
% snvQin:    SNV queue in 
%            1-d: top node in column 1, snv in column 2
%            2-d: top node row in column 1, column in column 2, snv in column 3
% top1,2:    top nodes of supernodes to be merged
%            NOTE:  assumes top1 is the PARENT sn of top2 sn
%            => in cssa alg, top2 has already been deleted from snvQ
% snvI:      SNV of the element to be inserted
%
% snout:     updated supernode data base
% snvQout:   updated SNV queue
%
% RGB INI September 1998
% RGB Rice January 1999

function [snout,snvQout] = condense(snin,snvQin,top1,top2);

sn = snin;

%----------------------------------------------------------------------------
if length(size(sn))==2       % 1-d case
  % values for update
  snvnew = ( sn(top1,3)*sn(top1,4) + sn(top2,3)*sn(top2,4) )/ ...
               ( sn(top2,4) + sn(top1,4) );
  numnew = sn(top1,4) + sn(top2,4);

  % point top node of the smaller sn to point to the top of the larger sn
  % also update data bases

  if sn(top1,4) >= sn(top2,4)    % top1 becomes the new top
    % point top2 to top1
    sn(top2,1) = top1;
  
    % update root pointers:   
    % root1 pointer keeps its -1 flag and top1 continues to point there 
    % (no changes)
    % root pointers at root2 and top2 set to zero
    if sn(top2,2) ~= -1          % => root2 and top2 are different nodes
    
      sn(sn(top2,2),2) = 0;
    end

    sn(top2,2) = 0; 
      
    % update snv
    sn(top1,3) = snvnew;
    % update num nodes  
    sn(top1,4) = numnew;
    % update snv queue by replacing top1 values with new top1 values
    % remember: top2 has already been deleted from snvQ
     snvQout = Qswap(snvQin,top1,top1,snvnew);
  
  else                          % top2 becomes the new top
    % point top1 to top2
    sn(top1,1) = top2;
 
    % update root pointers:   
    % root1 keeps its -1 root flag and top2 root pointer now points there
    if sn(top2,2) ~= -1           % => root2 and top2 are different nodes
      sn(sn(top2,2),2) = 0;
    end
    if sn(top1,2) ~= -1           % => root1 and top1 are different nodes
      sn(top2,2) = sn(top1,2);
      sn(top1,2) = 0;             % top1 now an internal node
    else                          % => root1 and top1 are the same node
    sn(top2,2) = top1;
    end  
  
    % update snv
    sn(top2,3) = snvnew;
    % update num nodes  
    sn(top2,4) = numnew;
    % update snv queue by replacing top1 values with new top2 values
    % remember: top2 has already been deleted from snvQ
    snvQout = Qswap(snvQin,top1,top2,snvnew);
   
  end
  
%----------------------------------------------------------------------------
else        % 2-d case
  top1R = top1(1); top1C = top1(2);
  top2R = top2(1); top2C = top2(2);
  
  % values for update
  snvnew = ( sn(top1R,top1C,5)*sn(top1R,top1C,6) + ...
             sn(top2R,top2C,5)*sn(top2R,top2C,6) )/ ...
                  ( sn(top1R,top1C,6) + sn(top2R,top2C,6) );
  numnew = sn(top1R,top1C,6) + sn(top2R,top2C,6);

  % point top node of the smaller sn to point to the top of the larger sn
  % also update data bases
  
  if sn(top1R,top1C,6) >= sn(top2R,top2C,6)    % top1 becomes the new top
    % point top2 to top1
    sn(top2R,top2C,1) = top1R;
    sn(top2R,top2C,2) = top1C;

    % update root pointers:   
    % root1 pointer keeps its -1 flag and top1 continues to point there 
    % (no changes)
    % root pointers at root2 and top2 set to zero
    if sn(top2R,top2C,3) ~= -1       % => root2 and top2 are different nodes
      sn( sn(top2R,top2C,3), sn(top2R,top2C,4), 3) = 0;
      sn( sn(top2R,top2C,3), sn(top2R,top2C,4), 4) = 0;
    end
    sn(top2R,top2C,3) = 0; 
    sn(top2R,top2C,4) = 0; 
      
    % update snv
    sn(top1R,top1C,5) = snvnew;
    % update num nodes  
    sn(top1R,top1C,6) = numnew;
    % update snv queue by replacing top1 values with new top1 values
    % remember: top2 has already been deleted from snvQ
    snvQout = Qswap(snvQin,top1,top1,snvnew);
  
  else                          % top2 becomes the new top
    % point top1 to top2
    sn(top1R,top1C,1) = top2R;
    sn(top1R,top1C,2) = top2C;
 
    % update root pointers:   
    % root1 keeps its -1 root flag and top2 root pointer now points there
    if sn(top2R,top2C,3) ~= -1    % => root2 and top2 are different nodes
      sn( sn(top2R,top2C,3), sn(top2R,top2C,4), 3) = 0;
      sn( sn(top2R,top2C,3), sn(top2R,top2C,4), 4) = 0;
    end
    
    if sn(top1R,top1C,3) ~= -1    % => root1 and top1 are different nodes
      sn(top2R,top2C,3) = sn(top1R,top1C,3);
      sn(top2R,top2C,4) = sn(top1R,top1C,4);
      sn(top1R,top1C,3) = 0;      % top1 now an internal node
      sn(top1R,top1C,4) = 0;      % top1 now an internal node
    else                          % => root1 and top1 are the same node
      sn(top2R,top2C,3) = top1R;
      sn(top2R,top2C,4) = top1C;
    end  
  
    % update snv
    sn(top2R,top2C,5) = snvnew;
    % update num nodes  
    sn(top2R,top2C,6) = numnew;
    % update snv queue by replacing top1 values with new top2 values
    % remember: top2 has already been deleted from snvQ
    snvQout = Qswap(snvQin,top1,top2,snvnew);
  end
end
snout = sn;
return



