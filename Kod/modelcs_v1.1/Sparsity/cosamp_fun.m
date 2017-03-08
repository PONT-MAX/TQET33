% cosamp_fun.m
% Performs CS reconstruction using the CoSaMP algorithm
% (D. Needell and J. Tropp, "CoSaMP: iterative signal recovery from 
%  incomplete and inaccurate measurements" , ACHA, 2009)
%
% Similar algorithm as cosamp.m . Uses function handles to implement
% multiplication with the matrix Phi
% 
% INPUTS
% yy           : measurements (M x 1)
% Phi_f, PhiT_f: function handles (refer cosamp_example for details)
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
% email:   chinmay@rice.edu

function [xhat,xcosamp] = cosamp_fun(yy, Phi_f, PhiT_f, N, K, Its);

%---
yy = yy(:); % 
M  = length(yy);
aa= zeros(N,Its); % stores current sparse estimate
s_cosamp = zeros(N,1);
kk=1; % current MP iteration
maxiter= 100;
verbose= 0;
tol= 1e-3;

while le(kk,Its),
    rr = yy - Phi_f(s_cosamp);
    proxy = PhiT_f(rr);
    %---Estimate support
    [tmp,ww]= sort(abs(proxy),'descend');
    tt= union(find(ne(s_cosamp,0)),ww(1:(2*K)));
    
    % Preparation for cg_solve
    PP_tt = @(z) A_I(Phi_f,z,tt,N);
    PP_transpose_tt = @(z) A_I_transpose(PhiT_f,z,tt);
    qq = PP_transpose_tt(yy);
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    %Pseudo-inverse
    [w, res, iter] = cgsolve(PPtranspose_PP_tt, qq, tol, maxiter, verbose);

    bb= 0*s_cosamp; bb(tt)= w;
    %---Prune
    kk = kk+1;
    [tmp,ww2]= sort(abs(bb),'descend'); 
    s_cosamp=0*bb;
    s_cosamp(ww2(1:K)) = bb(ww2(1:K));
    aa(:,kk) = s_cosamp;  
    
    if kk>1
       if norm(aa(:,kk)-aa(:,kk-1))<1e-2*norm(aa(:,kk))
           break;
       end
    end
end
xcosamp= aa;
xcosamp(:,kk:end)=[];
xhat = xcosamp(:,end);