function [ y_final, y_m ] = signal_process_48( filename, sidelength, ratio, f_s, show_plot  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5, show_plot = false; end

f_sw = 59.99*24 % Hz
N_yw = f_s/f_sw
M = round(sidelength^2*ratio/48)*48

% Read Data from file
y_struct = TDMS_getStruct(filename);
y_raw = y_struct.Measured_Data.AI0.data;
%figure(32); plot(y_raw)
% figure(99); plot(y_raw); pause();
% Automatic signal extraction
max_noise = abs(max(y_raw(1:f_s/10))*4) % From Pilot zero signal w/ buffer
mead_val = median(y_raw) % Will be Phi Signal stregth

% Find first Sensing matrix
i = 1;
while(true)
    if y_raw(i) > max_noise
        not_found = false;
        break;
    end
    i = i+1;
end


% When fist 100% pilot start is found all the first sensing matrix will
% allways be found 6 frames = N_sw*6*3
start_index = round(i + N_yw*6*24)

% find last sensing matrix
not_found = true;
i = start_index + 2000;

while(not_found)
    if y_raw(i) < mead_val*65/100
        not_found = false;
        break;
    end
    i = i+1;
    
end
while(y_raw(i) < min(y_raw(start_index:round(start_index+length(y_raw)*7/10))))
    i = i-1;
end

stop_index = i;
% 
 if show_plot
    figure(1); 
    plot(1:length(y_raw),y_raw,'b',[start_index stop_index],y_raw([start_index stop_index]),'*r');
 end

y = y_raw(start_index:stop_index);
y = y-min(y);


N_y = length(y);
N_yw = N_y/(M)
sm_window = zeros(1,N_y);
sm_window(1:N_yw:N_y) = max(y);
%
% Implement a lowpass filter

[b,a] = butter(2,0.3)
y_dn = filtfilt(b,a,y);

if true%show_plot
    figure(13)
    Y = fft(y);
    Y_dn = fft(y_dn);
    f = f_s*(0:(N_y-1))/N_y;
    subplot(2,1,1); plot(f,abs(Y)); axis([0 35000, 0 10000]); 
    title('Y');xlabel('f (Hz)');ylabel('|Y_gt(f)|')
    subplot(2,1,2); plot(f,abs(Y_dn)); axis([0 35000, 0 10000]); title('Ydn');xlabel('f (Hz)');ylabel('|Y_gt(f)|')
end



y_dn = medfilt1(y_dn,7,'truncate');

if true
    figure(22)
    subplot(2,1,2); plot(1:N_y,y_dn,'r',1:N_y,sm_window,'m',1:N_y,y,'b');axis([10000 11000, 0 max(y_dn)]);
    subplot(2,1,1);plot(1:N_y,y,'g',1:N_y,sm_window,'m');axis([10000 11000, 0 max(y)]);
end


NNN = round(N_yw/2.8)
% Gora denna kanslg for outliers
y_final = zeros(M,1);
yi = 1;
for i = 1:N_yw:N_y
    h = round(i+N_yw)-1;%round(i*N_yw)
    l = round(i);%ceil(h - N_yw)
    mf = mean(y_dn(l+round(N_yw/2.8):h-9));
    y_m(l:h) = mf;
    y_final(yi) = mf;
    yi = yi+1;
end

if show_plot
    figure(4)
    plot(1:N_y,y_dn,'g',1:N_y,y_m,'b')
end

end

