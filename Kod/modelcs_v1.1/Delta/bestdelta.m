% bestdelta.m
%
% Computes the best (K,delta,r)-approximation of a 1D signal
% Basically solves a linear program
%
% Needs CVX to run! To get CVX, go to:
% http://www.stanford.edu/~boyd/cvx/

function s = bestdelta(x,K,delta,r)

x=x(:); x = x';
N = length(x);
% set up linear programming equations
% objective cost
c = abs(x);
% all ones vector
a = ones(1,N);
% Windowing matrix
W = zeros(N-delta+1,N);
W(1,1:delta)=1;
for ii=2:N-delta+1
    W(ii,ii:ii+delta-1)=1;
end
cvx_begin
cvx_quiet(true)
    variable s(N)
    minimize(-c*s); 
    subject to 
        a*s <= K;
        W*s <= r;
        s >= 0;
        s <= 1;
cvx_end
