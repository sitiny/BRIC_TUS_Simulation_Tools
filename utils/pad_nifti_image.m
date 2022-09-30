function pad_nifti_image(input, output, pad_val)
%PAD_NIFTI_IMAGE Pads image to make space for transducer positioning
% 
% Inputs:
%   input:      Full filename and path to image you want to pad
%   output:     Full filename and path to where to save your output image
%   pad_val:    Integer or 1 x 3 array for number of voxels to pad by.
%               Padding will symmetrical on each side i.e. total padding =
%               pad_val * 2. If an integer is provided, padding will be the
%               same value in each dimension. 
% 
% Dependencies:
%   MATLAB Image Processing Toolbox for niftiinfo, niftiread, niftiwrite
% 
% Author: Siti N. Yaakub, University of Plymouth, 30 Sep 2022

% load image
img_info = niftiinfo(input);
img_vol = niftiread(input);

if length(pad_val) == 1
    pad_val = [pad_val, pad_val, pad_val];
    disp(['Padding in X Y Z by ' num2str(pad_val)])
elseif size(pad_val,1) == 1 && size(pad_val,2) == 3
    disp(['Padding in X Y Z by ' num2str(pad_val)])
elseif size(pad_val,1) == 3 && size(pad_val,2) == 1
    pad_val = pad_val';
    disp(['Padding in X Y Z by ' num2str(pad_val)])
else
    error('Pad value not valid.');
end
    
% pad images by pad_val on each side
new_vol = zeros(size(img_vol, 1) + pad_val(1) * 2, ...
    size(img_vol, 2) + pad_val(2) * 2, ...
    size(img_vol, 3) + pad_val(3) * 2);
new_vol = cast(new_vol, "like", img_vol);
new_vol(pad_val(1) + 1:size(img_vol, 1) + pad_val(1), ...
    pad_val(2) + 1:size(img_vol, 2) + pad_val(2), ...
    pad_val(3) + 1:size(img_vol, 3) + pad_val(3)) = img_vol;

% write new T1w MRI
img_info.Filename = []; img_info.Filemoddate = []; img_info.Filesize = [];
img_info.ImageSize = size(new_vol); 
tmp = img_info.Transform.T(4,:);
tmp = tmp - [pad_val(1)+1, pad_val(2)+1, pad_val(3)+1, 0];
img_info.Transform.T(4,:) = tmp;
niftiwrite(new_vol, output, img_info);