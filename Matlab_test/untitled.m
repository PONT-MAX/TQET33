for i = 0:7
2^i
a1 = zeros([256 256],'uint8');
a2 = zeros([256 256],'uint8');
a3 = zeros([256 256],'uint8');

a1 = a1 + (randi(2,[256 256],'uint8') - 1)*2^i;
a2 = a2 + (randi(2,[256 256],'uint8') - 1)*2^i;
a3 = a2 + (randi(2,[256 256],'uint8') - 1)*2^i;


end
%%
im_col = zeros([256 256,3],'uint8');
im_col(:,:,1) = a1;
im_col(:,:,2) = a2;
im_col(:,:,3) = a3;
figure(1); imshow(im_col(:,:,1), [0 255])

imwrite(im_col(:,:,1),'bw.png')