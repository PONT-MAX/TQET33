% Construct shapes
h = 2480;
w = 3511;

%% Cross
c = zeros(100,100,'uint8');
c(48:52,:) = 255; c(:,48:52) = 255;
c(46:54,[1:40,60:end]) = 255; c([1:40,60:end],46:54) = 255;
c(44:56,[1:30,70:end]) = 255; c([1:30,70:end],44:56) = 255;
c(42:58,[1:20,80:end]) = 255; c([1:20,80:end],42:58) = 255;
c(40:60,[1:10,90:end]) = 255; c([1:10,90:end],40:60) = 255;
figure(3); imshow(c, [0 255])
c2 = imresize(c,2,'bicubic');
figure(4); imshow(c2, [0 255])
%% Soft chess board
im = zeros(h,w,'uint8');

for x = 1:w
    for y = 1:h
        im(y,x) = 255*((2 + cos(x/80 + 3.14*2) + sin(y/80 + 3.14/2))/4);
    end
end

r_11 = 201-50; r_12 = r_11 + 199;
r_21 = 2211-50; r_22 = r_21 + 199;
c_11 = 206 -50; c_12 = c_11 + 199;
c_21 = 3218 - 50; c_22 = c_21 + 199;

im(r_11:r_12,c_11:c_12) = im(r_11:r_12,c_11:c_12) + c2;
im(r_11:r_12,c_21:c_22) = im(r_11:r_12,c_21:c_22) + c2;
im(r_21:r_22,c_11:c_12) = im(r_21:r_22,c_11:c_12) + c2;
im(r_21:r_22,c_21:c_22) = im(r_21:r_22,c_21:c_22) + c2;

r11 = 403:602; r22 = 1913:2112; c11 = 914:1113; c22 = 2413:2612;

im(r11,c11) = im(r11,c11) - c2;
im(r11,c22) = im(r11,c22) - c2;
im(r22,c11) = im(r22,c11) - c2;
im(r22,c22) = im(r22,c22) - c2;

figure(2); imshow(im, [0 255])
imwrite(im,'45.PNG','PNG');
%% Image 1, 10 shades of gray "Gray scale target (10):"

shades = 8;
im = zeros(h,w,'uint8');
step = 255/(shades-1);
ind_step = round(w/shades);
for i = 0:shades-1
    col = 255 - round(step*i);
    start_ind = i*ind_step + 1;
    stop_ind = start_ind + ind_step;
    im(:, start_ind:stop_ind ) = col; 
end
figure(2); imshow(im, [0 255])
imwrite(im,'1.PNG','PNG');


%% Image 3 Linear gradient

im = zeros(h,w,'uint8');
for i = 1:w
    col = 255 - floor(i/w*255);
    im(:, i ) = col; 
end
figure(2); imshow(im, [0 255])
imwrite(im,'2.PNG','PNG');

%% tilted rectangle

im = ones(h,w,'uint8')*217;

rec = zeros(h,h,'uint8');
rec(round(h/4):round(3*h/4),round(h/4):round(3*h/4)) = 153;
rec = imrotate(rec,30,'crop');
im(:,round(w/2-h/2):round(w/2+h/2)-1) = im(:,round(w/2-h/2):round(w/2+h/2)-1) - rec;
figure(1); imshow(im, [0 255])
imwrite(im,'3.PNG','PNG');

%% tilted rectangle inv

im = ones(h,w,'uint8')*64;

rec = zeros(h,h,'uint8');
rec(round(h/4):round(3*h/4),round(h/4):round(3*h/4)) = 153;
rec = imrotate(rec,30,'crop');
im(:,round(w/2-h/2):round(w/2+h/2)-1) = im(:,round(w/2-h/2):round(w/2+h/2)-1) + rec;
figure(1); imshow(im, [0 255])
imwrite(im,'4.PNG','PNG');

%% Triangles

im = ones(h, w, 'uint8')*217;
b = round(w/4);
t1 = zeros(b,b,'uint8');
yh = round((b-8)/2*cos(pi/6));
xCoords = [4 round((b-8)/2) b-4];
yCoords = [round((b-8)/2+yh) round((b-8)/2-yh) round((b-8)/2+yh)];
mask = poly2mask(xCoords, yCoords, b, b);
t1(mask) = 64;
t2 = imrotate(t1,90,'crop');
t3 = imrotate(t1,180,'crop');
t4 = imrotate(t1,270,'crop');
im(round(h/2-b/2):round(h/2+b/2)-1,:) = [t1 t2 t3 t4];
im(im == 0) = 217;

imwrite(im,'5.PNG','PNG');
figure(1); imshow(im, [0 255])

%% Inverted

im = ones(h, w, 'uint8')*64;
b = round(w/4);
t1 = zeros(b,b,'uint8');
yh = round((b-8)/2*cos(pi/6));
xCoords = [4 round((b-8)/2) b-4];
yCoords = [round((b-8)/2+yh) round((b-8)/2-yh) round((b-8)/2+yh)];
mask = poly2mask(xCoords, yCoords, b, b);
t1(mask) = 217;
t2 = imrotate(t1,90,'crop');
t3 = imrotate(t1,180,'crop');
t4 = imrotate(t1,270,'crop');
im(round(h/2-b/2):round(h/2+b/2)-1,:) = [t1 t2 t3 t4];
im(im == 0) = 64;

imwrite(im,'6.PNG','PNG');
figure(1); imshow(im, [0 255])

%% Squares

Im = zeros(h,w,'uint8');
paint = [1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1];
wy = round(w/5.1);
b = round(wy/12)
for row = 1:3
    for col = 1:5
        if paint(row,col) == 1
           Im(round(wy/3)+wy*(row-1)+b:round(wy/3)+wy*(row)-1,wy*(col-1)+b:wy*col-1 ) = 255; 
        end        
    end
end
figure(2); imshow(Im, [0 255])
imwrite(Im,'7.PNG','PNG');

%%
Im = zeros(h,w,'uint8');
paint = [1 1 0 0 1; 0 1 1 1 1; 1 1 0 0 1];
for row = 1:3
    for col = 1:5
        if paint(row,col) == 1
           Im(round(wy/3)+wy*(row-1)+b:round(wy/3)+wy*(row)-1,wy*(col-1)+b:wy*col-1 ) = 255; 
        end        
    end
end
figure(2); imshow(Im, [0 255])
imwrite(Im,'8.PNG','PNG');

%%
Im = zeros(h,w,'uint8');
paint = [0 0 1 0 0; 1 1 1 1 1; 0 0 1 0 0];
for row = 1:3
    for col = 1:5
        if paint(row,col) == 1
           Im(round(wy/3)+wy*(row-1)+b:round(wy/3)+wy*(row)-1,wy*(col-1)+b:wy*col-1 ) = 255; 
        end        
    end
end
figure(2); imshow(Im, [0 255])
imwrite(Im,'9.PNG','PNG');
%%
Im = zeros(h,w,'uint8');
paint = [1 1 1 1 1; 1 0 0 0 1; 1 1 1 1 1];
for row = 1:3
    for col = 1:5
        if paint(row,col) == 1
           Im(round(wy/3)+wy*(row-1)+b:round(wy/3)+wy*(row)-1,wy*(col-1)+b:wy*col-1 ) = 255; 
        end        
    end
end
figure(2); imshow(Im, [0 255])
imwrite(Im,'10.PNG','PNG');

%%
Im = zeros(h,w,'uint8');
paint = [0 0 0 0 0; 0 1 1 1 0; 0 0 0 0 0];
for row = 1:3
    for col = 1:5
        if paint(row,col) == 1
           Im(round(wy/3)+wy*(row-1)+b:round(wy/3)+wy*(row)-1,wy*(col-1)+b:wy*col-1 ) = 255; 
        end        
    end
end
figure(2); imshow(Im, [0 255])
imwrite(Im,'11.PNG','PNG');


%% BG

im = zeros(h,w,'uint8');
im(1:480,1:480) = 255;
im(end-480:end,1:480) = 255;
im(1:480,end-480:end) = 255;
im(end-480:end,end-480:end) = 255;
imwrite(im,'12.PNG','PNG');

