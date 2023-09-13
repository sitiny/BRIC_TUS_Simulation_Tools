function vox_coords = mm_to_vox(nifti_img, mm_coords)
% takes a 3x1 matrix of mm coordinates and converts to voxel coordinates
% requires nifti image file for affine transform
% returns a 3x1 matrix of voxel coordinates (round values as needed)

header=niftiinfo(nifti_img);
affine=header.Transform.T';
% check mmcoords is 3x1 rather than 1x3
if size(mm_coords,2) == 3
    mm_coords = mm_coords';
end

vox_coords = affine\[mm_coords;1];
vox_coords = vox_coords(1:3);
