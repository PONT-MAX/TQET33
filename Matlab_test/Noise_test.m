filename = '22_512.tdms';
sidelength = 512;
ratio = 0.3;
f_s = 69120;
N =  sidelength^2;
M = round(N*ratio/48)*48;
[y_dn, y_m] = signal_process_48( filename, sidelength, ratio, f_s, true  );
%%

%%

%dt_l = y_dn - y_final;
dt_mm = movmean(y_dn,50);
y_mm = y_dn - dt_mm;
y_mm = y_mm - min(y_mm);
figure(5)
plot(x,y_dn,'r',x,y_final,'b',x,dt_l,'g',x,y_mm,'g')

figure(6)
plot(x,dt_mm,'b',x,dt_l,'g')

%%
y_struct = TDMS_getStruct(filename);
y_raw = y_struct.Measured_Data.AI0.data;
figure(32); plot(y_raw)
%%
y_n = y_raw(1:120000);
figure(77);plot(y_n)
%%
m_n = mean(y_n)
std_n = std(y_n)
v_n = var(y_n)
y_nn = y_n - mean(y_n);
m_nn = mean(y_nn)
std_nn = std(y_nn)
v_nn = var(y_nn)
figure(71);plot(y_nn)
%%
y_s = y_m(1:120000) - mean(y_m(1:120000));
figure(78); plot(y_s)
dt_mm = movmean(y_s,50*48);
y_ss = y_s - dt_mm;
y_ss = y_ss - min(y_ss);
figure(79); plot(y_ss)
%%

xx = 1:120000;
figure(80);plot(xx,y_ss,'b',xx,(y_ss+y_nn-mean(y_nn)),'r')
xlabel('sample'); ylabel('Voltage [V]')
legend('Signal y', 'Background noise')
set(gca,'FontSize', 16)
%%
figure(81);plot(xx,y_ss,'b',xx,y_nn,'g')
xlabel('sample'); ylabel('Voltage [V]')
legend('Signal y', 'Background noise')
set(gca,'FontSize', 16)
%%
R = snr((y_ss-min(y_ss)),(y_nn-mean(y_nn)))
%%
figure(2); ylabel('y'); xlabel('sample')
set(gca,'fontsize',16)
figure(4); ylabel('y'); xlabel('sample')
set(gca,'fontsize',16)


%%
y_N = y_ss/max(y_ss);
var(y_N)
y_nnN = y_nn/max(y_ss);
var(y_nnN)
figure(21); plot(xx,y_N,'b',xx,y_nnN,'g')

