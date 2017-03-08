% deltarec.m
% Performs CS reconstruction of (K,\Delta)-signals
%

% %---IMPORTANT --- %
% This function calls bestdelta.m, which needs CVX to run. CVX is a Matlab 
% package for convex programming and is super easy to use. See: 
% http://www.stanford.edu/~boyd/cvx/
% for installation instructions.
% %--/Important --- %

% INPUTS
% yy      : measurements (M x 1)
% Phi     : measurement matrix (M x N)
% delta   : Minimum spike separation
% K       : Number of nonzeros
% Its     : number of iterations
% 
% OUTPUTS
% xhat    : Signal estimate (N x 1) 
% xdmp    : Matrix with N rows and at most Its columns; 
%           columns represent intermediate signal estimates   
% 
%
% CITE    : Chinmay Hegde, Marco Duarte, Volkan Cevher
%          "Compressive sensing of spike trains using a structured sparsity
%           model",SPARS, 2009.
% Created : Aug 2009.
% email   : chinmay@rice.edu

function [xhat,xdmp] = deltarec(yy, Phi, delta, K, Its);

yy = yy(:); % 
[M,N] = size(Phi);

xdmp = zeros(N,Its); 
xcosamp = zeros(N,Its);
kk=1; 
maxiter= 1000;
verbose= 0;
tol= 1e-3;
s_dmp = zeros(N,1);
while le(kk,Its),
    rdmp = yy - Phi*(s_dmp);
    proxy_dmp = Phi'*(rdmp);
    %--------------------------------------------%
    %-----Insert approximation algorithm here----%
    %--------------------------------------------%
    %-----(K,delta)-MP
    % Call bestdelta.m
    
    supp = bestdelta(proxy_dmp,2*K,delta,2);
    % check for weird cases
    flg = sum((supp<=.999).*(supp>=.001));
    if flg > 0
        disp('WARNING: nonbinary solutions to bestdelta');
    end
    newsupp = round(supp);
    tt_dmp=union(find(ne(s_dmp,0)),find(newsupp==1));  
    %-----/Insert----%

    %------Estimate------%
    [w_dmp, res, iter] = cgsolve(Phi(:,tt_dmp)'*Phi(:,tt_dmp), Phi(:,tt_dmp)'*yy,...
                                        tol,maxiter, verbose);
    bb1= zeros(N,1);
    bb1(tt_dmp)= w_dmp;

    %---Prune----%
    kk = kk+1;   
    
    %---Best Delta pruning of bb1
    supp = bestdelta(bb1,K,delta,1);
    % check for weird cases
    flg = sum((supp<=0.999).*(supp>=0.001));
    if flg > 0
        disp('WARNING: nonbinary solutions to bestdelta');
    end
    newsupp = round(supp);
    s_dmp=newsupp.*bb1;
    xdmp(:,kk) = s_dmp(:);

    if (norm(xdmp(:,kk)-xdmp(:,kk-1)) < 0.01*norm(xdmp(:,kk)))
        break;
    end
end
xdmp(:,kk+1:end)=[];
xhat = xdmp(:,end);