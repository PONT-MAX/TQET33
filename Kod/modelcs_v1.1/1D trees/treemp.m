% treemp.m
% Performs CS reconstruction of 1D wavelet-sparse signals using 
% the tree model
%
% IMPORTANT: Make sure you have the Rice Wavelet Toolbox in your path
%
% INPUTS
% yy      : measurements (M x 1)
% Phi     : measurement matrix (M x N)
% wvopts  : parameters for computing wavelet decomposition
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


function [xhat,xtmp] = treemp(yy, Phi, wvopts, K, Its);

yy = yy(:); % 
M  = length(yy);
N = wvopts.length;
L = wvopts.lv;
h = wvopts.filt;
%---
aa = zeros(N,Its); % stores current sparse estimate
kk=1; % current MP iteration
maxiter= 1000;
verbose= 0;
tol= 1e-3;
ws_tmp = zeros(N,1);

while le(kk,Its),
    rr = yy - Phi*(ws_tmp);
    proxy = Phi'*(rr);
    %-----Insert approximation algorithm here----%
    Bn = abs(proxy).^2;
    maskn = cssa1(Bn,L,2*K);  % computes best 2K term tree approximation
    % keep whichever value of maskn is nonzero
    tt=union(find(ne(ws_tmp,0)),find(maskn>0)); 
%     length(tt), pause(0.5)
    %-----/Insert----%

    
    %------Estimate------%
    [w, res, iter] = cgsolve(Phi(:,tt)'*Phi(:,tt), Phi(:,tt)'*yy,...
                                        tol,maxiter, verbose);
    bb= 0*ws_tmp; bb(tt)= w;
    
    %---Prune----%
    kk = kk+1;   
    % Do another CSSA approximation to retain best tree approximation
    Bn = abs(bb).^2;
    maskn = cssa1(Bn,L,K); 
    ws_tmp=0*bb; f=find(maskn>0);  ws_tmp(f) = bb(f);
    %[kk,length(f)]
    aa(:,kk) = midwt(ws_tmp,h,L);
    if kk>1
       if norm(aa(:,kk)-aa(:,kk-1))<5e-3*norm(aa(:,kk))
           break;
       end
    end
end
xtmp = aa;
xtmp(:,kk+1:end)=[];
xhat = xtmp(:,end);