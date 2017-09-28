function [ Im0, Im ] = get_swir_im( name, size )
%Get SWIR images, make sure that the image folder is in path

Im = importdata([int2str(name) '.tif']);
Im(:,513:end) = [];
Im = double(Im./(max(Im(:)/255)));
Im = imadjust(Im./255)*255;

if size == 128
Im0 = imresize(Im,0.25,'cubic');
Im = imresize(Im,0.5,'cubic');
Im(1:64,:) = [];
Im(:,1:64) = [];
Im(129:end,:) = [];
Im(:,129:end) = [];
elseif size == 256
Im0 = imresize(Im,0.5,'cubic');
Im(1:128,:) = [];
Im(:,1:128) = [];
Im(257:end,:) = [];
Im(:,257:end) = [];
elseif size == 512
    Im0 = Im;
end



end

