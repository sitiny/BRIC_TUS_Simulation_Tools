function CTX500_3D_skull_kArray_gen(varargin)
% =========================================================================
% CTX500_3D_skull_kArray_PCC_MRSTUS.m
%
% Run 3D k-wave acoustic simulations for the CTX-500 4-element transducer
% transducer modelled with kArray
% skull estimated from pCT generated from T1w-MRI
% target: PCC
%
% Siti N. Yaakub | 18-Aug-2022 17:11:03
%
% NB 30-Aug-2022 
% - fixed  idx2(:,2) = size(t1_img) + [shift_x, shift_x, shift_z]; --> 
% ...[shift_x, shift_y, shift_z];
% 
% - added option for (not) running thermal simulation
% - added option to run C++ as system command or from matlab
% - summary values are now written in a .csv file
% - images are displyed with the correct orientation
% =========================================================================

% do you want to run the thermal simulation? 
run_thermal_sim = false; 

% how to run the C++ simulation code
% http://www.k-wave.org/documentation/example_cpp_running_simulations.php
run_c = 'linux_system'; % possible values: 'matlab', 'terminal', 'linux_system';

% pathname for the input and output files
base_dir = '/mnt/nas8T/Experiments/TUS_social_comparison/01_Setup/TUS_simulation';
output_dir = fullfile(base_dir, 'Simulation');
load('driving_params_Isppa20.mat', 'driving_params');

if numel(varargin) == 0
    % subject-specific parameters
    subj_id         = 'golt0459';           % subject ID

    % add 1 if coordinates were obtained using fsleyes 
    % fsleyes uses zero-indexing but Matlab indexing starts at 1 rather than 0.
    focus_coords    = [91, 169, 174] + 1;   % point on the beam axis of the transducer [grid points]
    bowl_coords     = [101, 213, 211] + 1;  % centre of rear surface of transducer [grid points]
    focus_depth     = 46;                   % to nearest mm
else
    subj_id         = varargin{1};
    focus_coords    = varargin{2};
    bowl_coords     = varargin{3};
    focus_depth     = varargin{4};
end

% path to k-wave toolbox
% http://www.k-wave.org/installation.php
% http://www.k-wave.org/downloads/kWaveArray_alpha_0.3.zip
kWave_path = '/mnt/exp/utils/k-Wave';
addpath(genpath(kWave_path));

% --------------------
% SIMULATION VARIABLES
% --------------------

% medium parameters
c_min               = 1500;     % sound speed [m/s]
c_max               = 3100;     % max. speed of sound in skull (F. A. Duck, 2013.) [m/s]
rho_min             = 1000;     % density [kg/m^3]
rho_max             = 1900;     % max. skull density [kg/m3]
alpha_power         = 1.43;     % Robertson et al., PMB 2017 usually between 1 and 3? from Treeby paper
alpha_coeff_water   = 0;        % [dB/(MHz^y cm)] 0.05 from Fomenko et al., 2020?
alpha_coeff_min     = 4;        % 
alpha_coeff_max     = 8.7;      % [dB/(MHz cm)] Fry 1978 at 0.5MHz: 1 Np/cm (8.7 dB/cm) for both diploe and outer tables

% source parameters
freq        = 500e3;            % source frequency [Hz]
source_roc  = 63.2e-3;          % bowl radius of curvature [m]

% pulse parameters
pulse_length    = 20e-3;	% pulse length [s]
pulse_rep_freq  = 5;        % pulse repetition frequency [Hz]
stim_dur        = 80;       % duration of stimulation [s]

% aperture diameters of the elements given an inner, outer pairs [m]
% 0.0254 = conversion from inches to m (element dimensions given in inches)
diameters       = [0 1.28; 1.3 1.802; 1.82 2.19; 2.208 2.52].' * 0.0254;

% grid parameters
hu_min 	= 300;	% minimum Hounsfield Unit in CT image
hu_max 	= 2000;	% maximum Hounsfield Unit in CT image

% computational parameters
ppw             = 3;        % number of points per wavelength
record_periods  = 3;        % number of periods to record
cfl             = 0.3;      % CFL number

%% prepare simulation

% --------------------
% SKULL MODEL
% --------------------

% Load CT image (nifti format)
% voxel size = 1 x 1 x 1 mm3, matrix size: varies, mostly 176x256x256
t1_img = niftiread(fullfile(base_dir, 'PseudoCT', 'T1toN4', ['sub-' subj_id '_T1w_n4_masked.nii.gz']));
input_ct = niftiread(fullfile(base_dir, 'PseudoCT', 'N4topCT', ['sub-' subj_id '_T1w_n4_masked_pct.nii.gz']));
header = niftiinfo(fullfile(base_dir, 'PseudoCT', 'N4topCT', ['sub-' subj_id '_T1w_n4_masked_pct.nii.gz']));
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
idx2(:,2) = size(t1_img) + [shift_x, shift_y, shift_z];
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

%%
tic;
% driving parameters (depends on focus depth)
idx = find(driving_params.dist == focus_depth);
source_amp      = repmat(driving_params.amp{idx}, [1,4]);   % source pressure [Pa]
source_phase    = deg2rad([driving_params.phase{idx,:}]);   % phase [rad]

% create time varying source
source_sig = createCWSignals(kgrid.t_array, freq, source_amp, source_phase);

% assign source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

% --------------------
% SENSOR
% --------------------

% set sensor mask (restrict mask to head mask based on t1_img - assumes t1_img has been masked for pct creation)
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(t1_img>0) = 1;

% record the pressure
sensor.record = {'p'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * ppp + 1;

%%% Run simulation

% --------------------
% ACOUSTIC SIMULATION
% --------------------

% set input options
input_args = {'PMLSize', 10, 'PMLInside', true, 'PlotPML', true, 'DisplayMask', source.p_mask};

if strcmp(run_c, 'matlab') % run C++ in Matlab

    sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});

elseif strcmp(run_c, 'terminal') % run C++ in terminal

% % input and output filenames (these must have the .h5 extension)
% input_filename  = fullfile(output_dir, ['CERMEP_' subj_id '_CTX500_3D_35kPa_kArray_PCC_input.h5']);
% output_filename = fullfile(output_dir, ['CERMEP_' subj_id '_CTX500_3D_35kPa_kArray_PCC_output.h5']);
% sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', input_filename);
% % display the required syntax to run the C++ simulation
% disp(['Using a terminal window, navigate to the ' filesep 'binaries folder of the k-Wave Toolbox']);
% disp('Then, use the syntax shown below to run the simulation:');
% disp(['./kspaceFirstOrder-OMP -i ' input_filename ' -o ' output_filename ' --p_raw -s ' num2str(sensor.record_start_index)]);
% to re-input results into matlab
% sensor_data.p = h5read(output_filename, '/p');

elseif strcmp(run_c, 'linux_system') % run C++ using system commands from matlab
    % define names of temporary input and output files
    input_filename = fullfile(output_dir, ['kwave_input_' subj_id '.h5']);
    output_filename = fullfile(output_dir, ['kwave_output_' subj_id '.h5']);
    
    % create input file
    sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', input_filename);

    % run simulation using bash commands
    kspacefct_path = fullfile(kWave_path, 'binaries', 'kspaceFirstOrder-OMP');
    [status, simul_output] = system([kspacefct_path ' -i ' input_filename ...
       ' -o ' output_filename ' --p_raw -s ' num2str(sensor.record_start_index) '-t 12'], '-echo'); %#ok<ASGLU> 

    if status
        error('simulation failed');
    end

    % re-input results into matlab
    sensor_data.p = h5read(output_filename, '/p');
    
    % delete temporary input and output files
    delete input_filename output_filename
end
toc;

% --------------------
%% CALCULATE PRESSURE
% --------------------

% extract amplitude from the sensor data
p = extractAmpPhase(sensor_data.p, 1/kgrid.dt, freq, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
% p = reshape(p, Nx, Ny, Nz);
tmp_p = p;
p = zeros(size(model));
p(t1_img>0) = tmp_p;
clear tmp_p;

% calculate acoustic intensities
[max_pressure, idx] = max(p(:)); % [Pa]
[mx, my, mz] = ind2sub(size(p), idx);

if model(mx, my, mz) > 1
    Isppa = max_pressure^2 / (2 * max(medium.density(:)) * max(medium.sound_speed(:))); % [W/m2]
elseif model(mx, my, mz) == 0
    Isppa = max_pressure^2 / (2 * rho_min * c_min); % [W/m2]
end
Isppa = Isppa * 1e-4; % [W/cm2]
Ispta = Isppa * pulse_length * pulse_rep_freq; % [W/cm2]

% MI = max_pressure (in MPa) / sqrt freq (in MHz)
MI = max_pressure * 1e-6 / sqrt(freq * 1e-6);

% find -6dB focal volume
focal_vol = length(find(p > 0.5*max(p(:))));

% calculate FWHM in each axis
field_profile_x = p(:,my,mz);
field_profile_y = p(mx,:,mz);
field_profile_z = squeeze(p(mx,my,:));

plot_fwhm=0;
field_focal_size_6dB(1) = fwhm2(field_profile_x, dx, mx, plot_fwhm);
field_focal_size_6dB(2) = fwhm2(field_profile_y, dx, my, plot_fwhm);
field_focal_size_6dB(3)= fwhm2(field_profile_z, dx, mz, plot_fwhm);

p_focus = p(focus_coords(1), focus_coords(2), focus_coords(3));
isppa_focus = p_focus^2 / (2 * rho_min * c_min) * 1e-4;

% --------------------
%% SUMMARY
% --------------------
clear source;
disp(subj_id)
disp(['Shift = [' num2str(shift_x) ', ' num2str(shift_y) ', ' num2str(shift_z) ']'])
disp(['PPW = ' num2str(ppw)])
disp(['CFL = ' num2str(cfl)])
disp(['Coordinates of max pressure: [' num2str(mx-shift_x) ', ' num2str(my-shift_y) ', ' num2str(mz-shift_z) ']'])
disp(['# Voxels > 90% max pressure = ' num2str(length(find(p > (0.9*max_pressure))))])
disp(['# Voxels > 50% max pressure = ' num2str(length(find(p > (0.5*max_pressure))))])
disp(['Distance from transducer rear surface: ' num2str(norm(bowl_coords-[mx,my,mz])*dx*1e3) ' mm'])
disp(['Max Pressure = ' num2str(max_pressure * 1e-6) ' MPa'])
disp(['MI = ' num2str(MI)])
disp(['Isppa = ' num2str(Isppa) ' W/cm2'])
disp(['Ispta = ' num2str(Ispta*1e3) ' mW/cm2'])
disp(['Pressure at focus = ' num2str(p_focus * 1e-6) ' MPa']);
disp(['Isppa at focus = ' num2str(isppa_focus) ' W/cm2'])
disp(['-6dB focal volume = ' num2str(focal_vol) ' mm3'])
disp(['-6dB max width = ' num2str(max(field_focal_size_6dB)*1e3) ' mm'])
disp(' ')

save(fullfile(output_dir, [subj_id '_CTX500_3D_kArray_5Hz_rTUS.mat']));

%% copying results to csv file
% create result file if it does not exist and write header
if ~exist(fullfile(output_dir, 'simulation_results.csv'), 'file')
 disp('Result file does not exist, creating file.')
 fileID = fopen(fullfile(output_dir, 'simulation_results.csv'), 'w' );  
 fprintf(fileID, '%s\n', ...
    ['id,focus_coord, bowl_coord, focus_depth, PPW, CFL, Coordinates of max pressure,' ...
    '# Vxls > 90% max pressure', '# Vxls > 50% max pressure',...
    ' Distance from transducer rear surface (mm), Max Pressure (mPA), MI, Isppa (W/cm2),' ...
    ' Pressure at focus (MPa), Isppa at focus ( W/cm2), -6dB focal volume (mm3), -6dB max width (mm)']);
 fclose(fileID);
end 

% write values
fileID = fopen(fullfile(output_dir, 'simulation_results.csv'),'a');
fprintf(fileID,['%s, %d %d %d, %d %d %d, %d, ' ...
    '%f, %f, %d %d %d, %d, ' ...
    '%d, %f, ' ...
    '%f, %f, %f, %f, %f, %f, %f\n'], ...
    subj_id, focus_coords, bowl_coords, focus_depth,  ...
    ppw, cfl, mx-shift_x, my-shift_y, mz-shift_z, length(find(p > (0.9*max_pressure))), ...
    length(find(p > (0.5*max_pressure))), norm(bowl_coords-[mx,my,mz])*dx*1e3, ...
    max_pressure * 1e-6, MI, Isppa, p_focus * 1e-6, isppa_focus, focal_vol, max(field_focal_size_6dB)*1e3);
fclose(fileID);

% --------------------
%% VISUALISATION
% --------------------

figure;
ax1a = axes;
imagesc(ax1a, imrotate(squeeze(t1_img(mx,:,:)), 90), [50,500]);
% axis square; 
hold all;
ax2a = axes;
im2 = imagesc(ax2a, imrotate(squeeze(p(mx,:,:))*1e-6, 90));
im2.AlphaData = 0.5;
% axis square;
linkaxes([ax1a,ax2a]); ax2a.Visible = 'off'; ax2a.XTick = []; ax2a.YTick = []; 
colormap(ax1a,'gray')
colormap(ax2a,'turbo')
set([ax1a,ax2a],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2a,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1a,'Acoustic Pressure Amplitude overlaid on T1')
saveas(gcf, fullfile(output_dir, [subj_id '_5Hz_rTUS_sag.jpg']));

figure;
ax1b = axes;
imagesc(ax1b, imrotate(squeeze(t1_img(mx,:,:)), 90), [50,500]);
% axis square;
hold all;
ax2b = axes;
im2 = imagesc(ax2b, imrotate(squeeze(p(mx,:,:)>(0.5*max_pressure))*1e-6, 90));
im2.AlphaData = 0.5;
% axis square;
linkaxes([ax1b,ax2b]); ax2b.Visible = 'off'; ax2b.XTick = []; ax2b.YTick = []; 
colormap(ax1b,'gray')
colormap(ax2b,'turbo')
set([ax1b,ax2b],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2b,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1b,'50% Acoustic Pressure Amplitude')
saveas(gcf, fullfile(output_dir, [subj_id '_5Hz_rTUS_sag_50%.jpg']));

figure;
ax1c = axes;
imagesc(ax1c, imrotate(squeeze(t1_img(mx,:,:)), 90), [50,500]);
% axis square;
hold all;
ax2c = axes;
im2 = imagesc(ax2c, imrotate(squeeze(p(mx,:,:)>(0.9*max_pressure))*1e-6, 90));
im2.AlphaData = 0.5;
% axis square;
linkaxes([ax1c,ax2c]); ax2c.Visible = 'off'; ax2c.XTick = []; ax2c.YTick = []; 
colormap(ax1c,'gray')
colormap(ax2c,'turbo')
set([ax1c,ax2c],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2c,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1c,'90% Acoustic Pressure Amplitude')
saveas(gcf, fullfile(output_dir, [subj_id '_5Hz_rTUS_sag_90%.jpg']));

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(:,my,:)), 90), [50,500]);
% axis square; 
hold all;
ax2a = axes;
im2 = imagesc(ax2a, imrotate(squeeze(p(:,my,:))*1e-6, 90));
im2.AlphaData = 0.5;
% axis square;
linkaxes([ax1,ax2a]); ax2a.Visible = 'off'; ax2a.XTick = []; ax2a.YTick = []; 
colormap(ax1,'gray')
colormap(ax2a,'turbo')
set([ax1,ax2a],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2a,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude overlaid on T1')
saveas(gcf, fullfile(output_dir, [subj_id '_5Hz_rTUS_cor.jpg']));

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(:,:,mz)), 90), [50,500]);
% axis square; 
hold all;
ax2a = axes;
im2 = imagesc(ax2a, imrotate(squeeze(p(:,:,mz))*1e-6, 90));
im2.AlphaData = 0.5;
% axis square;
linkaxes([ax1,ax2a]); ax2a.Visible = 'off'; ax2a.XTick = []; ax2a.YTick = []; 
colormap(ax1,'gray')
colormap(ax2a,'turbo')
set([ax1,ax2a],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2a,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude overlaid on T1')
saveas(gcf, fullfile(output_dir, [subj_id '_5Hz_rTUS_ax.jpg']));

%
%% In case the max presure falls outside the brain (redo summary and plots if running this)
%
% find max p in brain zz from plot of pressure in sagittal plane - x-axis coord below current max P

zz = 154;
tmp=p(:,:,1:zz); 

[max_pressure, idx] = max(tmp(:)); % [Pa]
[mx, my, mz] = ind2sub(size(tmp), idx);

if model(mx, my, mz) > 1
    Isppa = max_pressure^2 / (2 * max(medium.density(:)) * max(medium.sound_speed(:))); % [W/m2]
elseif model(mx, my, mz) == 0
    Isppa = max_pressure^2 / (2 * rho_min * c_min); % [W/m2]
end
Isppa = Isppa * 1e-4; % [W/cm2]
Ispta = Isppa * pulse_length * pulse_rep_freq; % [W/cm2]

% MI = max_pressure (in MPa) / sqrt freq (in MHz)
MI = max_pressure * 1e-6 / sqrt(freq * 1e-6);

% find -6dB focal volume
focal_vol = length(find(tmp > 0.5*max(tmp(:))));

% calculate FWHM in each axis
field_profile_x = tmp(:,my,mz);
field_profile_y = tmp(mx,:,mz);
field_profile_z = squeeze(tmp(mx,my,:));

plot_fwhm=0;
field_focal_size_6dB(1) = fwhm2(field_profile_y, dx, my, plot_fwhm);
field_focal_size_6dB(2) = fwhm2(field_profile_x, dx, mx, plot_fwhm);
field_focal_size_6dB(3)= fwhm2(field_profile_z, dx, mz, plot_fwhm);

% % --------------------
%% THERMAL SIMULATION
% % --------------------

if ~run_thermal_sim
   return 
end

% convert the absorption coefficient to nepers/m
alpha_np = db2neper(medium.alpha_coeff, medium.alpha_power) * ...
    (2 * pi * freq).^medium.alpha_power;

% reshape the data, and calculate the volume rate of heat deposition
Q = alpha_np .* p.^2 ./ (medium.density .* medium.sound_speed);

run thermal simulation
tic;
% clear the input structures
% clear medium source sensor;
clear source;

% set the background temperature and heating term
source2.Q = Q;
source2.T0 = 37;

% define medium properties related to diffusion
% ref: https://itis.swiss/virtual-population/tissue-properties/database
% in skull
medium2.density = rho_min + (rho_max - rho_min) * (model - 0) / (hu_max - 0);
medium2.thermal_conductivity = zeros(size(model));
medium2.thermal_conductivity(model > 0) = 0.32;
medium2.specific_heat = 826 + (2524 - 826) * (1-((model - hu_min) / (hu_max - hu_min)));
% in water
medium2.density(model == 0)              = rho_min;     % [kg/m^3]
medium2.thermal_conductivity(model == 0) = 0.6;      % [W/(m.K)]
medium2.specific_heat(model == 0)        = 4178;     % [J/(kg.K)]

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium2, source2, [], 'DisplayUpdates', false, 'PlotSim', false);

% set source on time and off time
on_time  = pulse_length;  % [s]
off_time = 1/pulse_rep_freq - pulse_length;  % [s]
num_bursts = stim_dur * pulse_rep_freq;
% num_bursts = 2;
% macaque protocol
% pulse_length = 30e-3;
% pulse_rep_freq = 10;
% stim_dur = 40;
% on_time  = pulse_length;  % [s]
% off_time = 1/pulse_rep_freq - pulse_length;  % [s]
% num_bursts = stim_dur * pulse_rep_freq;

% set time step size
% dt = 0.1;
dt = on_time/1;
% dt = on_time/4;

maxT1 = zeros(size(model));
for nb = 1:num_bursts
    disp(['burst #' num2str(nb)])
    % turn on heat source
    kdiff.Q = Q;
    % take time steps
    kdiff.takeTimeStep(round(on_time / dt), dt);
    
    % store the current temperature field
    T1 = kdiff.T;
    
    % turn off heat source and take time steps
    kdiff.Q = 0;
    kdiff.takeTimeStep(round(off_time / dt), dt);
    
    % store the current temperature field
    T2 = kdiff.T;
end
toc;

disp(['Max temperature rise = ' num2str(max(T1(:))-source2.T0) ' degC'])
save(fullfile(output_dir, [subj_id '_CTX500_3D_kArray_5Hz_rTUS_therm.mat']));

% --------------------
%% VISUALISATION
% --------------------

% plot the thermal dose and lesion map (yz)
figure;

% plot the acoustic pressure
subplot(2, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.z_vec * 1e3, squeeze(p(mx,:,:) * 1e-6).');
h = colorbar;
xlabel(h, '[MPa]');
ylabel('z-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Acoustic Pressure Amplitude');

% plot the volume rate of heat deposition
subplot(2, 2, 2);
imagesc(kgrid.y_vec * 1e3, kgrid.z_vec * 1e3, squeeze(Q(mx,:,:) * 1e-7).');
h = colorbar;
xlabel(h, '[kW/cm^2]');
ylabel('z-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Volume Rate Of Heat Deposition');

% plot the temperature after heating
subplot(2, 2, 3);
imagesc(kgrid.y_vec * 1e3, kgrid.z_vec * 1e3, squeeze(T1(mx,:,:)).');%, [37, 37.1]);
h = colorbar;
xlabel(h, '[degC]');
ylabel('z-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Temperature After Heating');

% plot the temperature after cooling
subplot(2, 2, 4);
imagesc(kgrid.y_vec * 1e3, kgrid.z_vec * 1e3, squeeze(T2(mx,:,:)).');%, [37, 37.1]);
h = colorbar;
xlabel(h, '[degC]');
ylabel('z-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Temperature After Cooling');

% set colormap and enlarge figure window
% colormap(jet(256));
colormap('turbo');
scaleFig(2, 2);
saveas(gcf, fullfile(output_dir, [subj_id '_5Hz_rTUS_therm_yz.jpg']));

% plot the thermal dose and lesion map (xz)
figure;

% plot the acoustic pressure
subplot(2, 2, 1);
imagesc(kgrid.x_vec * 1e3, kgrid.z_vec * 1e3, squeeze(p(:,my,:) * 1e-6).');
h = colorbar;
xlabel(h, '[MPa]');
ylabel('z-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Acoustic Pressure Amplitude');

% plot the volume rate of heat deposition
subplot(2, 2, 2);
imagesc(kgrid.x_vec * 1e3, kgrid.z_vec * 1e3, squeeze(Q(:,my,:) * 1e-7).');
h = colorbar;
xlabel(h, '[kW/cm^2]');
ylabel('z-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Volume Rate Of Heat Deposition');

% plot the temperature after heating
subplot(2, 2, 3);
imagesc(kgrid.x_vec * 1e3, kgrid.z_vec * 1e3, squeeze(T1(:,my,:)).');%, [37, 37.1]);
h = colorbar;
xlabel(h, '[degC]');
ylabel('z-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Temperature After Heating');

% plot the temperature after cooling
subplot(2, 2, 4);
imagesc(kgrid.x_vec * 1e3, kgrid.z_vec * 1e3, squeeze(T2(:,my,:)).');%, [37, 37.1]);
h = colorbar;
xlabel(h, '[degC]');
ylabel('z-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Temperature After Cooling');

% set colormap and enlarge figure window
% colormap(jet(256));
colormap('turbo');
scaleFig(2, 2);
saveas(gcf, fullfile(output_dir, [subj_id '_5Hz_rTUS_therm_xz.jpg']));

% plot the thermal dose and lesion map (xy)
figure;

% plot the acoustic pressure
subplot(2, 2, 1);
imagesc(kgrid.x_vec * 1e3, kgrid.y_vec * 1e3, (p(:,:,mz) * 1e-6).');
h = colorbar;
xlabel(h, '[MPa]');
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Acoustic Pressure Amplitude');

% plot the volume rate of heat deposition
subplot(2, 2, 2);
imagesc(kgrid.x_vec * 1e3, kgrid.y_vec * 1e3, (Q(:,:,mz) * 1e-7).');
h = colorbar;
xlabel(h, '[kW/cm^2]');
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Volume Rate Of Heat Deposition');

% plot the temperature after heating
subplot(2, 2, 3);
imagesc(kgrid.x_vec * 1e3, kgrid.y_vec * 1e3, T1(:,:,mz).');%, [37, 37.1]);
h = colorbar;
xlabel(h, '[degC]');
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Temperature After Heating');

% plot the temperature after cooling
subplot(2, 2, 4);
imagesc(kgrid.x_vec * 1e3, kgrid.y_vec * 1e3, T2(:,:,mz).');%, [37, 37.1]);
h = colorbar;
xlabel(h, '[degC]');
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Temperature After Cooling');

% set colormap and enlarge figure window
% colormap(jet(256));
colormap('turbo');
scaleFig(2, 2);
saveas(gcf, fullfile(output_dir, [subj_id '_5Hz_rTUS_therm_xy.jpg']));
