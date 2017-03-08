% treemp_fun.m
%
% Similar algorithm as treemp.m . Uses function handles to implement
% multiplication with the matrix Phi
% 
% INPUTS
% yy           : measurements (M x 1)
% Phi_f, PhiT_f: function handles (refer cosamp_example for details)
% wvopts       : parameters for computing wavelet decomposition
% K            : signal sparsity
% N            : signal length
% Its          : number of iterations
% 
% OUTPUTS
% xhat   : Signal estimate (N x 1) 
% xcosamp: Matrix with N rows and at most Its columns; 
%          columns represent intermediate signal estimates   
% 
%
% CITE:    Richard Baraniuk, Volkan Cevher, Marco Duarte, Chinmay Hegde
%          "Model-based compressive sensing", submitted to IEEE IT, 2008.
% Created: Aug 2009.

function [xhat,xtmp] = treemp_fun(yy, Phi_f, PhiT_f, wvopts, K, Its);

%---
yy = yy(:); % 
M  = length(yy);
N = wvopts.length;
L = wvopts.lv;
h = wvopts.filt;

aa= zeros(N,Its); % stores current sparse estimate
ws_tmp = zeros(N,1);
kk=1; % current MP iteration
maxiter= 1000;
verbose= 0;
tol= 1e-3;

while le(kk,Its),
    rr = yy - Phi_f(ws_tmp);
    proxy = PhiT_f(rr);
    %---Estimate support
    Bn = abs(proxy).^2;
    maskn = cssa1(Bn,L,2*K); 
    tt=union(find(ne(ws_tmp,0)),find(maskn>0)); 
    %length(tt), pause(0.5)
    
    % Preparation for cg_solve
    PP_tt = @(z) A_I(Phi_f,z,tt,N);
    PP_transpose_tt = @(z) A_I_transpose(PhiT_f,z,tt);
    qq = PP_transpose_tt(yy);
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    %Pseudo-inverse
    [w, res, iter] = cgsolve(PPtranspose_PP_tt, qq, tol, maxiter, verbose);

    bb= 0*ws_tmp; bb(tt)= w;
    %---Prune
    kk = kk+1;
    Bn = abs(bb).^2;
    maskn = cssa1(Bn,L,K); 
    ws_tmp=0*bb; f=find(maskn>0);  ws_tmp(f) = bb(f);
    aa(:,kk) = midwt(ws_tmp,h,L);  
    %kk
    if kk>1
       if norm(aa(:,kk)-aa(:,kk-1))<5e-3*norm(aa(:,kk))
           break;
       end
    end
end
xtmp = aa;
xtmp(:,kk+1:end)=[];
xhat = xtmp(:,end);