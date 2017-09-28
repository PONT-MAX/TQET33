function [ Im ] = generate_pattern( size)
%Input is 64,128 or 256 and generates an test image w/ res 64x64 || 128x128
%|| 256x256 w/ some pattern on it

Im = zeros(256,256);
Im(10:100,10:100) = 255.0;
Im(110:100,110:180) = 255.0;
Im(200:250,15:65) = 255.0;
Im(150:200,150:220) = 255.0;
Im(20:120,200:230) = 255.0;

if size == 128
    scale = 0.5;
    Im = imresize(double(Im),scale);
elseif size == 64
    scale = 0.25;
    Im = imresize(double(Im),scale);
end
    

end

