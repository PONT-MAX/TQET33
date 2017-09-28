y = rand(1,10)/10 + 5;
ys = zeros(1,5000);
k = 1;
for i = 1:10
    for j = 1:500
       ys(k) = y(i) + 0.0035*rand(1);
       k = k + 1;
    end

end


figure(1)
plot(ys)

Ys1 = fft(ys);
figure(4)
plot(abs(Ys1))

x = 1:5000;
ys = ys + 0.003*cos(x/180);



figure(2)
plot(ys)

Ys2 = fft(ys);
figure(3)
plot(abs(Ys2))

figure(2)
plot(abs(Ys1-Ys2))