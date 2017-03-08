% CSSA-based tree matching pursuit : no functions!
% Chin: July 16, '08

function xhat= cosamp_2D_nf(yy, Phi, wvopts, K, Its, wn);

yy = yy(:); % 
M  = length(yy);
n = wvopts.xlength;
N = wvopts.slength;
L = wvopts.lv;
h = wvopts.filt;
%---
aa = zeros(N,Its); % stores current sparse estimate
kk=1; % current MP iteration
maxiter= 1000;
verbose= 0;
tol= 1e-3;
ws_vec = zeros(N,1);
ws_im = zeros(n,n);
while le(kk,Its),
    rr = yy - Phi*(ws_vec);
    proxy = Phi'*(rr);
    %proxy_im = reshape(proxy, n, n);
    %-----Insert approximation algorithm here----%

    [tmp,ww]= sort(abs(proxy),'descend');
    tt = union(find(ne(ws_vec,0)),ww(1:(2*K)));
    
    %------Estimate------%
    [w, res, iter] = cgsolve(Phi(:,tt)'*Phi(:,tt), Phi(:,tt)'*yy,...
                                        tol,maxiter, verbose);
    bb=zeros(N,1);
    bb(tt)= w;
    
    %---Prune----%
    kk = kk+1;   
    [tmp,ind]= sort(abs(bb),'descend');
    ws_vec = 0*bb;
    ws_vec(ind(1:K)) = bb(ind(1:K));
    ws_im = reshape(ws_vec, n, n);
%     figure(4), hold on, stem(ws_vec, 'k*', 'MarkerSize', 0.5)
%     stem(wn, 'rd', 'MarkerSize', 0.5), hold off, shg, pause(0.1)
    [kk, sum(abs(ws_vec)>0), norm(ws_vec-wn)]
    x_im = midwt(ws_im,h,L); aa(:,kk) = x_im(:);
    figure(3), imagesc(x_im), colormap(gray), axis image, rmaxis, shg, pause(.01)
%     hold on, plot(aa(:,kk), 'b'), hold off
    if kk>1
       if norm(aa(:,kk)-aa(:,kk-1))<1e-2*norm(aa(:,kk))
           break;
       end
    end
end
xhat = aa;
xhat(:,kk:end)=[];