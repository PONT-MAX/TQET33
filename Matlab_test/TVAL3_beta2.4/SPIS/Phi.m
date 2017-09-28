function [ A ] = Phi( N , ratio)
%Returns Sensing matrix N = 64, 128, 256, 512
% User have to extract M from A after ths function.

if N == 64^2 || N == 128^2
    A = double(fWHtrans(eye(N))*N);
    A = (A + 1)./2;
    
    % Do the row permutation direct (same for all)
    load(['row_perm_', int2str(sqrt(N)), '.mat'])
    load(['col_perm_', int2str(sqrt(N)), '.mat'])
    
    col_perm_wh = zeros(1,N);
    col_perm_wh(col_perm) = 1:N;
    
    if ratio > 1 || ratio < 0.05 % If stupid
        ratio = 0.2
    end
    
    picks = row_perm(1:N-1);
    M = round(N*ratio/3)*3;
    
    A = A(picks(1:M),:);
    A = A(:,col_perm_wh);

    
elseif N == 256^2
    A = single(Phi_256(ratio));
    
elseif N == 512^2
    A = Phi_512(ratio);
else
    A = zeros(1,1);
    'WRONG SIZE try w/ N = 64^2, 128^2, 256^2 or 512^2'
end

end

