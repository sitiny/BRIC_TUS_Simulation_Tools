%% prepare simulation

% --------------------
% SKULL MODEL
% --------------------

% Load CT image (nifti format)
% voxel size = 1 x 1 x 1 mm3, matrix size: varies, mostly 176x256x256
input_ct = niftiread(fullfile(base_dir, 'Sims', 'pct', ['sub-' subj_id '_T1w_n4_masked_pct.nii']));
t1_img = niftiread(fullfile(base_dir, 'Sims', 't1', ['sub-' subj_id '_T1w_n4_masked.nii']));
header = niftiinfo(fullfile(base_dir, 'Sims', 'pct', ['sub-' subj_id '_T1w_n4_masked_pct.nii']));
input_ct = double(input_ct);

% update hu_max
ct_max = max(input_ct(:));
if ct_max < hu_max
    hu_max = ct_max;
end
clear ct_max;

% truncate CT HU (see Marsac et al., 2017)
skull_model = input_ct;
skull_model(skull_model < hu_min) = 0; % only use HU for skull acoustic properties
skull_model(skull_model > hu_max) = hu_max;
% skull_model(focus_coords(1),focus_coords(2),focus_coords(3))=10000;

% pad images by 100 on each side
tmp_model = zeros(size(skull_model,1)+200, size(skull_model,2)+200, size(skull_model,3)+200);
tmp_model(101:size(skull_model,1)+100, ...
    101:size(skull_model,2)+100, ...
    101:size(skull_model,3)+100) = skull_model;
tmp_focus = focus_coords+100;

% centre on focus
% grid size = 256x256x256, new focus coords = [128,128,128]
shift_idx = [tmp_focus(1)-128+1,tmp_focus(1)+128;
    tmp_focus(2)-128+1,tmp_focus(2)+128; ...
    tmp_focus(3)-128+1,tmp_focus(3)+128];
model = tmp_model(shift_idx(1,1):shift_idx(1,2), ...
    shift_idx(2,1):shift_idx(2,2), ...
    shift_idx(3,1):shift_idx(3,2));

shift_x = 128-focus_coords(1);
shift_y = 128-focus_coords(2);
shift_z = 128-focus_coords(3);

bowl_coords     = bowl_coords + [shift_x, shift_y, shift_z];	% centre of rear surface of transducer [grid points]
focus_coords    = focus_coords + [shift_x, shift_y, shift_z];  % point on the beam axis of the transducer [grid points]

% move t1 image
new_t1 = zeros(size(model));
idx = zeros(3,2);
for ii = 1:3
    if shift_idx(ii,1)-100 < 0
        idx(ii,1) = 1;
    elseif shift_idx(ii,1)-100 > 0
        idx(ii,1) = shift_idx(ii,1)-100;
    end
    if shift_idx(ii,2)-100 <= size(t1_img,ii)
        idx(ii,2) = shift_idx(ii,2)-100;
    elseif shift_idx(ii,2)-100 > size(t1_img,ii)
        idx(ii,2) = size(t1_img,ii);
    end
end
idx2 = zeros(3,2);
idx2(:,1) = [1,1,1] + [shift_x, shift_y, shift_z];
idx2(:,2) = size(t1_img) + [shift_x, shift_x, shift_z];
for ii = 1:3
    if idx2(ii,1) < 0
        idx2(ii,1) = 1;
    end
    if idx2(ii,2) > size(model,ii)
        idx2(ii,2) = size(model,ii);
    end
end
tmp_t1 = t1_img(idx(1,1):idx(1,2), idx(2,1):idx(2,2), idx(3,1):idx(3,2));
new_t1(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2)) = ... 
    t1_img(idx(1,1):idx(1,2), idx(2,1):idx(2,2), idx(3,1):idx(3,2));
t1_img = new_t1;
clear new_t1;

% --------------------
% MEDIUM
% --------------------

% assign medium properties for skull
% derived from CT HU based on Marsac et al., 2017 & Bancel et al., 2021
medium.density = rho_min + (rho_max - rho_min) * ... 
                (model - 0) / (hu_max - 0);
medium.sound_speed = c_min + (c_max - c_min) * ... 
                    (medium.density - rho_min) / (rho_max - rho_min);
medium.alpha_coeff = alpha_coeff_min + (alpha_coeff_max - alpha_coeff_min) * ...
                    (1 - (model - hu_min) / (hu_max - hu_min)).^0.5;

% assign medium properties for non-skull (brain, soft tissue, modelled as water)
medium.density(model == 0) = rho_min;
medium.sound_speed(model == 0) = c_min;
medium.alpha_coeff(model == 0) = alpha_coeff_water;

medium.alpha_power = alpha_power;

% --------------------
% GRID
% --------------------

% calculate the grid spacing based on the PPW and F0
dx = c_min / (ppw * freq);   % [m]
% % calculate ppw based on grid/image spacing and F0
% ppw = c_min / (d * freq);

% compute the size of the grid
[Nx, Ny, Nz] = size(model);     % [grid points]

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% compute points per temporal period
ppp = round(ppw / cfl);

% compute corresponding time spacing
dt = 1 / (ppp * freq);
dt_stability_limit = checkStability(kgrid, medium);
if dt_stability_limit ~= Inf
    dt = dt_stability_limit;
end

% calculate the number of time steps to reach steady state
t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2) / c_min; 

% create the time array using an integer number of points per period
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
ppw = c_min / (dx * freq);
cfl = c_min * dt / dx;
disp(['PPW = ' num2str(ppw)]);
disp(['CFL = ' num2str(cfl)]);

% --------------------
% SOURCE
% --------------------

% set bowl position and orientation
bowl_pos = [kgrid.x_vec(bowl_coords(1)), ...
            kgrid.y_vec(bowl_coords(2)), ...
            kgrid.z_vec(bowl_coords(3))];
focus_pos = [kgrid.x_vec(focus_coords(1)), ...
             kgrid.y_vec(focus_coords(2)), ...
             kgrid.z_vec(focus_coords(3))];

% create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.01, 'UpsamplingRate', 16);

% add bowl shaped element
karray.addAnnularArray(bowl_pos, source_roc, diameters, focus_pos)

% assign binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% visualise transducer (check it is within the grid)
model_mask = model;
model_mask(model_mask>0)=1;
tmp = source.p_mask+model_mask;
volumeViewer(tmp);
clear model_mask tmp;
% figure; 
% subplot(1,3,1); imagesc(squeeze(model(focus_coords(1),:,:)));
% subplot(1,3,2); imagesc(squeeze(model(:,focus_coords(2),:)));
% subplot(1,3,3); imagesc(squeeze(model(:,:,focus_coords(3))));