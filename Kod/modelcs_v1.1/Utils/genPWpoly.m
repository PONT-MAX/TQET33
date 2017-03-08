function y = genPWpoly(N,maxOrder,breakpoints);
% N - signal length
% maxOrder - degree of polynomials
% breakpoints - number of breakpoints (located at random)

b = 1+randperm(N-1);
b = sort(b(1:breakpoints));

coeffs = randn(maxOrder+1,breakpoints+1);

leftVec = [1 b];
rightVec = [b-1 N];

y = zeros(N,1);

for segment = 1:breakpoints+1  
  left = leftVec(segment);
  right = rightVec(segment);
  
  mid = mean([left right]);
  
  xVals = ([left:right]-mid)*2/N;
  
  yVals = zeros(size(xVals));
  for order = 0:maxOrder
    yVals = yVals + coeffs(order+1,segment)*xVals.^order;
  end
  
  y(left:right) = yVals;
  
end