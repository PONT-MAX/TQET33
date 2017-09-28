% Record or test
clear all; close all
record = true
% Choose Image size and scaling
sidelength = 256;
N = sidelength^2;
scale = 512/sidelength;

% Ratio
ratio = 0.20;
M = round(N*ratio/3)*3; % 3 / frame
%


% Get perm  Walsh Hadamard matrix
WH = Phi(N,ratio);
WH(1:10,1:10)

%
if record
    v = VideoWriter([int2str(sidelength),'_M_', int2str(M),'_r_', int2str(ratio*100) ,'.avi'],'Uncompressed AVI');
    open(v)
end

% Init current frame to write
Im = zeros(sidelength*scale,sidelength*scale,3,'uint8');

% start w/ 3 frames blac, 3 frames 100%, 2 framers black
for i = 1:8
    if i == 4
        Im = Im + 255;
    end
    if i == 7
        Im = Im*0;
    end
    %['Frame: ', int2str(i), ', Maxval: ', int2str(max(Im(:)))]
    if record
        writeVideo(v,Im);
        ['Frame: ', int2str(i), ' written, FC = ', int2str(v.FrameCount)]
        figure(2);
        subplot(1,3,1); imshow(Im(:,:,1)); title(['i: ',int2str(i)])
        subplot(1,3,2); imshow(Im(:,:,2))
        subplot(1,3,3); imshow(Im(:,:,3))
        figure(3); imshow(Im)
        pause(1)
    end
end

% Before reshape matrix (holder)
Im_br = zeros(1,N,1,'uint8');
% Before resize matrix (holder)
Im_bs = zeros(sidelength,sidelength,1,'uint8');

%First Sensing matrix sohult apper in two frames due to rise time from zero
%Encode into last ch, will make signal consistent
Im_br = uint8(WH(1,:));
Im_bs = reshape(Im_br,sidelength,sidelength);
Im(:,:,3) = imresize(Im_bs,scale, 'nearest')*255;
Im(:,:,1) = Im(:,:,1)*0;
Im(:,:,2) = Im(:,:,2)*0;
% Write frame

if record
    writeVideo(v,Im);
    ['Frame: ', int2str(9), ' written, FC = ', int2str(v.FrameCount)]
    figure(2);
    subplot(1,3,1); imshow(Im(:,:,1)); title(['i: ',int2str(9)])
    subplot(1,3,2); imshow(Im(:,:,2))
    subplot(1,3,3); imshow(Im(:,:,3))
    figure(3); imshow(Im)
    pause(1)
end

for frame = 1:3:M
    if mod(frame, 100) == 0 % Progress write out
        frame
    end
    % Ch 1
    Im_br = uint8(WH(frame,:));
    Im_bs = reshape(Im_br,sidelength,sidelength);
    Im(:,:,1) = imresize(Im_bs,scale, 'nearest')*255;
    
    % Ch 2
    Im_br = uint8(WH(frame + 1,:));
    Im_bs = reshape(Im_br,sidelength,sidelength);
    Im(:,:,2) = imresize(Im_bs,scale, 'nearest')*255;
    
    % Ch 3
    Im_br = uint8(WH(frame + 2,:));
    Im_bs = reshape(Im_br,sidelength,sidelength);
    Im(:,:,3) = imresize(Im_bs,scale, 'nearest')*255;
    
    % Write frame
    if record
        writeVideo(v,Im);
    end
    
end
'After WH encode'
fc = v.FrameCount
fc*3 - 9*3 - M

% Add 10 empty frames last
Im = Im*0;
for i = 1:10
    if record
        writeVideo(v,Im);
    end
end

if record
    close(v)
end

