% treemp_nf_2D.m
% Performs CS reconstruction of 2D wavelet-sparse signals using 
% the quad-tree model
%
% IMPORTANT: Make sure you have the Rice Wavelet Toolbox in your path
% http://dsp.rice.edu/software/rice-wavelet-toolbox
%
% differs from treemp_greedy_2D in choice of tree-approximation alg
%
% INPUTS
% yy      : measurement (M x 1)
% Phi     : measurement matrix (M x N)
% wvopts  : parameters for computing wavelet decomposition
%           check treemp_2D_example for details
% K       : signal sparsity
% Its     : number of iterations
% 
% OUTPUTS
% xhat    : Signal estimate (N x 1) 
% xtmp :    Matrix with N rows and at most Its columns; 
%           columns represent intermediate signal estimates   
% 
%
% CITE    : Richard Baraniuk, Volkan Cevher, Marco Duarte, Chinmay Hegde
%          "Model-based compressive sensing", submitted to IEEE IT, 2008.
% Created : Aug 2009.
% email   : chinmay@rice.edu

function xhat= treemp_nf_2D(yy, Phi, wvopts, K, Its, wn);

path(path,'../Utils')
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
    proxy_im = reshape(proxy, n, n);
    %-----Insert approximation algorithm here----%

    Bn = abs(proxy_im).^2;
%     dbstop if error
    if (norm(size(Bn)-[n,n])>0)
        break;
    end
    maskn = cssa2(Bn,L,2*K); 
    % keep whichever value of maskn is nonzero
    tt=union(find(ne(ws_vec,0)),find(maskn(:)>0)); 
%     length(tt), pause(0.5)
    %-----/Insert----%

    %------Estimate------%
    [w, res, iter] = cgsolve(Phi(:,tt)'*Phi(:,tt), Phi(:,tt)'*yy,...
                                        tol,maxiter, verbose);
    bb=zeros(N,1);
    bb(tt)= w;
    
    %---Prune----%
    kk = kk+1;   
    % Do another CSSA approximation to retain best tree approximation
    bb_im = reshape(bb,n,n);
    Bn = abs(bb_im).^2;
    disp('Pruning')
    maskn = cssa2(Bn,L,K); 
    ws_im=0*bb_im; f=find(maskn>0); 
    ws_im(f)=bb_im(f); ws_vec=ws_im(:);
    [kk,length(f)]
    x_im = midwt(ws_im,h,L); aa(:,kk) = x_im(:);
    figure(1), imagesc(x_im), colormap(gray), axis image, shg, pause(.1)
%     hold on, plot(aa(:,kk), 'b'), hold off
    if kk>1
       if norm(aa(:,kk)-aa(:,kk-1))<1e-2*norm(aa(:,kk))
           break;
       end
    end
end
xhat = aa;
xhat(:,kk:end)=[];