tmp_focusbin = p > 0.5*max(p(:));
tmp_focusbin = int16(tmp_focusbin);

idx = zeros(3,2);
idx(1,:) = [1, size(skull_model,1)]; % no shift laterally (x, dim=1), just centre focus
for dim = 2:3 % shift to make space anterior (y, dim=2) and superior (z, dim=3)
    if shiftdim(dim,2) < size(skull_model,dim)
        idx(dim,:) = [1, shiftdim(dim,2)];
    else
        idx(dim,:) = [shiftdim(dim,1),size(skull_model,dim)];
    end
end

focal_vol_bin = zeros(size(input_ct),'int16');
focal_vol_bin(idx(1,1):idx(1,2), idx(2,1):idx(2,2), idx(3,1):idx(3,2)) = ...
    tmp_focusbin(padidx(1,1):padidx(1,2), padidx(2,1):padidx(2,2), padidx(3,1):padidx(3,2));

header.Filename=[]; header.Filemoddate=[]; header.Filesize=[]; 
header.Datatype='int16'; header.BitsPerPixel=16;
niftiwrite(focal_vol_bin, fullfile(output_dir, ['sub-' subj_id '_6dBfocalvolume.nii']), header);

%%
tmp_focusbin = p > 0.5*max(p(:));
tmp_focusbin = int16(tmp_focusbin);

idx = zeros(3,2);
for ii = 1:3
    if shift_idx(ii,1)-100 < 0
        idx(ii,1) = 1;
    elseif shift_idx(ii,1)-100 > 0
        idx(ii,1) = shift_idx(ii,1)-100;
    end
    if shift_idx(ii,2)-100 <= size(input_ct,ii)
        idx(ii,2) = shift_idx(ii,2)-100;
    elseif shift_idx(ii,2)-100 > size(input_ct,ii)
        idx(ii,2) = size(input_ct,ii);
    end
end
idx2 = zeros(3,2);
idx2(:,1) = [1,1,1] + [shift_x, shift_y, shift_z];
idx2(:,2) = size(input_ct) + [shift_x, shift_x, shift_z];
for ii = 1:3
    if idx2(ii,1) < 0
        idx2(ii,1) = 1;
    end
    if idx2(ii,2) > size(model,ii)
        idx2(ii,2) = size(model,ii);
    end
end

focal_vol_bin = zeros(size(input_ct),'int16');
focal_vol_bin(idx(1,1):idx(1,2), idx(2,1):idx(2,2), idx(3,1):idx(3,2)) = ...
    tmp_focusbin(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2));

header.Filename=[]; header.Filemoddate=[]; header.Filesize=[]; 
header.Datatype='int16'; header.BitsPerPixel=16;
niftiwrite(focal_vol_bin, fullfile(output_dir, ['sub-' subj_id '_6dBfocalvolume.nii']), header);

% 
% tmp_t1 = t1_img(idx(1,1):idx(1,2), idx(2,1):idx(2,2), idx(3,1):idx(3,2));
% new_t1(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2)) = ... 
%     t1_img(idx(1,1):idx(1,2), idx(2,1):idx(2,2), idx(3,1):idx(3,2));