
k = -20*pi:0.1:20*pi;
a = atan(k);
a = a-min(a);
a = a/max(a);
n = length(a);
N = (1:n)./(max(n)*0.125);
i = ones(1,n);
i(1,1:end/2) = 0;
n90 = ones(1,n)*0.9;
n10 = ones(1,n)*0.1;
figure(1); plot(N,a,'b',N,i,'m',N,n90,'r',N,n10,'r')
axis([1, 7, -0.1 1.1]); xlabel('pixels'); ylabel('intensity')
legend('cameras measured edge', 'ideal edge', '10% and 90% intensity bounds')
set(gca,'FontSize', 16);