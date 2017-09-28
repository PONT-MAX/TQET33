
samp = [20:27 30:37 42:43 60:65 70 76 79];
n = length(samp);
mean_volt = zeros(1,n);
sigma = zeros(1,n);
snr_v = zeros(1,n);
k = 1;
for i = samp
nr = int2str(i)
filename = ['signal/', nr, '_512.tdms'];
sidelength = 512;
ratio = 0.3;
f_s = 69120;
N =  sidelength^2;
M = round(N*ratio/48)*48;
[y_dn, y_m] = signal_process_48( filename, sidelength, ratio, f_s, false  );
figure(99); plot(y_dn)
%dt_l = y_dn - y_final;
dt_mm = movmean(y_dn,50);
y_mm = y_dn - dt_mm;
y_mm = y_mm - min(y_mm);
%figure(5)
%plot(x,y_dn,'r',x,y_final,'b',x,dt_l,'g',x,y_mm,'g')

%figure(6)
%plot(x,dt_mm,'b',x,dt_l,'g')


y_struct = TDMS_getStruct(filename);
y_raw = y_struct.Measured_Data.AI0.data;
%figure(32); plot(y_raw)
if k == 11
xx = 1:12000;
else
xx = 1:120000;
end
y_n = y_raw(xx);
%figure(77);plot(y_n)

m_n = mean(y_n);
std_n = std(y_n);
v_n = var(y_n);
y_nn = y_n - mean(y_n);
m_nn = mean(y_nn);
std_nn = std(y_nn);
v_nn = var(y_nn);
%figure(71);plot(y_nn)
%
y_s = y_m(xx) - mean(y_m(xx));
%figure(78); plot(y_s)
dt_mm = movmean(y_s,50*48);
y_ss = y_s - dt_mm;
y_ss = y_ss - min(y_ss);
%figure(79); plot(y_ss)
%



figure(80);plot(xx,y_ss,'b',xx,(y_ss+y_nn-mean(y_nn)),'r')
xlabel('sample'); ylabel('Voltage [V]')
legend('Signal y', 'Background noise')
set(gca,'FontSize', 16)
%

figure(81);plot(xx,y_ss,'b',xx,y_nn,'g')
xlabel('sample'); ylabel('Voltage [V]')
legend('Signal y', 'Background noise')
set(gca,'FontSize', 16)
%

SNR_1 = snr((y_ss-min(y_ss)),(y_nn-mean(y_nn)))

%figure(2); ylabel('y'); xlabel('sample')
%set(gca,'fontsize',16)
%figure(4); ylabel('y'); xlabel('sample')
%set(gca,'fontsize',16)


y_N = y_ss/max(y_ss);
var(y_N);
y_nnN = y_nn/max(y_ss);
variance_noise = var(y_nnN)
mean_voltage = median(y_raw)

figure(21); plot(xx,y_N,'b',xx,y_nnN,'g')

mean_volt(k) = mean_voltage;
sigma(k) = variance_noise;
snr_v(k) = SNR_1;
k = k + 1;
end

%%
good = [3 4 5 10 11 21 22 23 24 25];
half = [6 14 15 16 17 18 19 20];
bad = [1 2 7 8 9 12 13 26 27];

s_g = snr_v(good)
s_h = snr_v(half)
s_b = snr_v(bad)

mv_g = mean_volt(good)
mv_h = mean_volt(half)
mv_b = mean_volt(bad)

figure(1);
for k = good
   plot(mv_g,s_g,'*g','MarkerSize', 14); hold on
end

for k = half
   plot(mv_h,s_h,'*y','MarkerSize', 14); hold on
end

for k = bad
   plot(mv_b,s_b,'*r','MarkerSize', 14) ; hold on
end


%plot(mean_volt,sigma,'x','MarkerSize', 10)
%figure(2); plot(snr_v,mean_volt,'.')
xlabel('Mean signal intensity [V]'); ylabel('SNR')
set(gca,'fontsize',16,'LineWidth',0.5)


