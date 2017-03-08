
% function maskp(x,w,nlevels,winsize,dB);
%
% Inputs:
%   x - 1-D signal
%   w - wavelet transform of x
%   nlevels - number of levels in wavelet transform
%   winsize - data window display size 
%     (plot displays x(1:winsize), 
%      to hide wraparound effect) 
%   dB - number of dB for image plot     * BUT WE DO LINEAR SCALE
%
% Outputs:
%   Plot of signal and energy in its TF tiling
%
% Subroutines:
%   setind.m for indexing wavelet coefficients
%
% NOTE: This is basically Matt's WAVEPLOT.m (my DWTPLOT.m)
% modified by
% me and renamed because of the PC's restriction to
% 7 character file names (!)

function maskp(x,w,nlevels,winsize,dB);

N=length(w);

% Sizing vector w to a row vector
w=reshape(w,1,length(w));

% TF tiling of wavelet coefficients
TFout=zeros(2^nlevels,N);

lind=setind(nlevels,N);

% Form dc tiles 
dcind=1:(min(lind)-1);
wband=repmat(w(dcind), round(2^nlevels),1);
wband=reshape(wband,1,N);
TFout(1,:)=wband;

% Forming ac tiles
for i=nlevels:-1:1
  lind=setind(i,N);
  time_inc = 2^i;
  f_inc = round(2^(nlevels-i));
  f = round(1+2^(nlevels-i));
  wband=repmat(w(lind),time_inc,1);
  wband=reshape(wband,1,N);
  wband=repmat(wband,f_inc,1);
  TFout(f:f+f_inc-1,:)=wband;
end


% Can use other colormaps here
%colormap('gray')

% Can use other monotonic functions of nrg here
% nrgout=log10(1.0+TFout.^2);
nrgout=log10(1.0+abs(TFout));

% TAKE ABS VALUE OF TFOUT!

tfimage(abs(flipud(TFout)),'li',dB) 

return







subplot('position',[0 .8 1 .2]);
plot(x);
axis([0 winsize-1 min(x) max(x)])
axis('off')
subplot('position',[0 0 1 .79]);
imagesc(-nrgout(:,1:winsize))
axis('off')


%set(imhan,'axesXTickMode','off','axesYTickMode','off');
%xhan=xlabel('time','FontSize',18)
%yhan=ylabel('frequency','FontSize',18);













