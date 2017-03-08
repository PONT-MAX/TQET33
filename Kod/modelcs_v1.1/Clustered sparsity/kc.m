% kc.m
% Performs CS reconstruction of 1D clustered-sparse signals
%

% %---IMPORTANT --- %
% This function makes use of code written in C++, which you'll need to mex.
% I've included pre-compiled binaries for Mac OS X and 32-bit Linux;
% for other platforms, do as follows:
% 1) navigate to the 'Utils' folder in the Matlab command window
% 2) type: 'mex bestKCsparse_l2mex.cpp'
% 3) If no error shows up, you're good to go.
% 
% %--/Important --- %

% INPUTS
% yy      : measurements (M x 1)
% Phi     : measurement matrix (M x N)
% K       : Number of nonzeros
% C       : number of clusters
% Its     : number of iterations
% 
% OUTPUTS
% xhat    : Signal estimate (N x 1) 
% xkc     : Matrix with N rows and at most Its columns; 
%           columns represent intermediate signal estimates   
% 
%
% CITE    : Volkan Cevher, Piotr Indyk, Chinmay Hegde, Richard Baraniuk
%          "Recovery of clustered-sparse signals from compressive
%           measurements",SampTA, 2009.
% Created : Aug 2009.
% email   : chinmay@rice.edu
%
% Hat tip : Drew Waters for the C++ code.


function [xhat,xkc]= kc(yy, Phi, K, C, Its);

path(path,'../Utils')

[M,N] = size(Phi);
aa = zeros(N,Its);
kk=1;
maxiter= 1000;
verbose= 0;
tol= 1e-3;
tic
while  le(kk,Its),
    %--- Calculate the residual
    rr= yy- Phi*aa(:,kk);
    proxy = Phi'*rr;
   
    %---best (2K,2C)-estimate of the residual ---%
    tt= bestKCsparse_l2mex(abs(proxy),2*K,2*C);
    tt= union(find(ne(aa(:,kk),0)),find(tt>0));
    
    %--- signal estimation
    bb= zeros(N,1);
    bb(tt)= pinv(Phi(:,tt))*yy;

    %--best (K,C)-estimate of the estimate 
    tt2= bestKCsparse_l2mex(abs(bb),K,C);
    tt2= find(tt2>0);
    kk= kk+1;
    aa(tt2,kk)= bb(tt2);
    %---
    [kk-1,toc]
    if (norm(aa(:,kk)-aa(:,kk-1)) < 1e-2*norm(aa(:,kk)))
       break;
    end    
end
xkc = aa;
xkc(:,kk+1:end)=[];
xhat = xkc(:,end);