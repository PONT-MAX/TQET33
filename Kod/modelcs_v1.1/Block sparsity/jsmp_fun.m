% jsmp_fun.m
%
% Similar algorithm as jsmp.m . Uses function handles to implement
% multiplication with the matrix Phi
% 
% INPUTS
% yy           : measurements (M x 1)
% Phi_f, PhiT_f: function handles (refer jsmp_example for details)
% J            : block length
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
%
% The model for block-sparsity used here is equivalent to the JSM-2 
% model introduced in:
% Baron et al, "Distributed Compressive Sensing", preprint 2009


function [xhat,xjsmp] = jsmp_fun(yy, Phi_f, PhiT_f, N, K, J, Its);

%---
yy = yy(:); % 
M  = length(yy);
aa= zeros(N,Its); % stores current sparse estimate
s_jsmp = zeros(N,1);
num_blocks = round(N/J);

kk=1; % current MP iteration
maxiter= 100;
verbose= 0;
tol= 1e-3;

while le(kk,Its),
    rr = yy - Phi_f(s_jsmp);
    proxy = PhiT_f(rr);
    %---Estimate support
   % pick the K largest blocks. Easy.
    proxy_jsmp_block = reshape(proxy, J, num_blocks);
    [trash,blockww] = sort(sum(proxy_jsmp_block.^2,1),'descend');
    newsupp = zeros(J,num_blocks);
    newsupp(:,blockww(1:(2*K))) = 1;
    newsupp = reshape(newsupp, N, 1);
    tt=union(find(ne(s_jsmp,0)),find(newsupp==1));  
    
    % Preparation for cg_solve
    PP_tt = @(z) A_I(Phi_f,z,tt,N);
    PP_transpose_tt = @(z) A_I_transpose(PhiT_f,z,tt);
    qq = PP_transpose_tt(yy);
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    %Pseudo-inverse
    [w, res, iter] = cgsolve(PPtranspose_PP_tt, qq, tol, maxiter, verbose);

    bb1= 0*s_jsmp; bb1(tt)= w;
    %---Prune
    kk = kk+1;
    bb1_block = reshape(bb1, J, num_blocks);
    [trash,blockww] = sort(sum(bb1_block.^2,1),'descend');
    newsupp = zeros(J,num_blocks);
    newsupp(:,blockww(1:K)) = 1;
    newsupp = reshape(newsupp, N, 1);
    s_jsmp=0*s_jsmp;
    s_jsmp = bb1.*newsupp;
    aa(:,kk) = s_jsmp;
    
    if kk>1
       if norm(aa(:,kk)-aa(:,kk-1))<5e-3*norm(aa(:,kk))
           break;
       end
    end
end
xjsmp= aa;
xjsmp(:,kk+1:end)=[];
xhat = xjsmp(:,end);