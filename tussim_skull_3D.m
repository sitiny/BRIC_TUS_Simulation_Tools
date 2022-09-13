function tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, ...
    focus_coords, bowl_coords, focus_depth, transducer, varargin)
%TUSSIM_SKULL_3D Run 3D k-wave acoustic simulation transcranially with skull
%   estimated using CT (or pseudo-CT) images.
%
% The simulation grid is 256x256x256, with voxel size of 1x1x1mm^3. The
% driving parameters are only valid for the CTX-500 with free-field Isppa
% of 20 W/cm^2. If you would like the simulation for a different free-field
% Isppa or a different transducer, you will have to provide your own source
% pressure and phase as optional Name-Value paired input arguments.
%
% You will need to supply a co-registered T1-weighted MR image and CT (or
% pseudo-CT) image for use in the simulations. These images must have voxel
% sizes of 1x1x1mm^3. The matrix size can be any size and will be
% automatically adjusted in the script to allow space for placing your
% transducer within the grid.
%
% It is preferable that you provide a T1-weighted MRI that has had noise
% outside the head masked out (e.g. the one used for pseudo-CT generation)
% so that the sensor is set to within this head mask. Alternatively, you
% can provide a brain extracted T1-weighted MRI instead, but this means you
% will not be able to simulate temperature rise at the skull interface when
% running the thermal simulation.
%
% Running the script without the acoustic or thermal simulation allows you
% to check that the transducer position relative to the head is correct.
% If you are satisfied, re-run the script with the option:
% 'RunAcousticSim', true.
%
% Currently, there is no functionality to run only a thermal simulation
% without running the acoustic simulation. Alternatively you can load the
% saved acoustic simulation .mat file and run the thermal simulation cells.
%
% Usage:
%   tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, ...
%       focus_coords, bowl_coords, focus_depth, transducer)
%   tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, ...
%       focus_coords, bowl_coords, focus_depth, transducer, Name, Value)
%   tussim_skull_3D('sub-test01', 'sub-test01_t1w.nii', ...
%       'sub-test01_pct.nii', 'output_dir', ...
%       [99, 161, 202], [90, 193, 262], 60, ...
%       'CTX500', 'RunAcousticSim', true);
%
% Inputs:
%   subj_id:        ID of the subject you are running the simulation for.
%   t1_filename:    Full file path to the T1-weighted MR image.
%   ct_filename:    Full file path to the CT (or pseudo-CT) image.
%   output_dir:     Full path to the output directory.
%   focus_coords:   3-element array of voxel coordinates of the desired TUS
%                   focus. Add 1 if reading these off a viewer that uses
%                   zero-indexing (MATLAB indexing starts from 1).
%   bowl_coords:    3-element array of voxel coordinates of the centre of
%                   the transducer base. Add 1 if reading these off a
%                   viewer that uses zero-indexing.
%   focus_depth:    Distance from transducer face to intended focus in mm,
%                   rounded to the nearest integer.
%   transducer:     Options are 'CTX500', 'CTX560' or 'CTX250'.
%
% Optional Name-Value pair inputs:
%   'RunAcousticSim':   Boolean controlling whether to run acoustic
%                       simulation (default: false).
%   'RunThermalSim':    Boolean controlling whether to run thermal
%                       simulation(default: false).
%   'PulseLength':      Pulse length in s (default: 20e-3 s).
%   'PulseRepFreq':     Pulse repetition frequency in Hz (default: 5 Hz).
%   'StimDuration':     Duration of TUS in s (default: 80 s).
%   'SourcePressure':   Source pressure in Pa.
%   'SourcePhase':      4-element array of phases of each transducer
%                       element in degrees for the focal depth required.
%   'RunCppCode':       Where to run the C++ simulation code (default: 'matlab')
%                       Options are 'matlab', 'terminal', 'linux_system'
%
% WARNING: Default acoustic values in function are for 500 kHz transducers.
% Please set your own values if using a transducer with a central frequency
% other than 500kHz.
%
% Dependencies:
%   k-Wave Toolbox (http://www.k-wave.org)
%   kArray (http://www.k-wave.org/downloads/kWaveArray_alpha_0.3.zip)
%     	Copyright (C) 2009-2017 Bradley Treeby
%
% Author: Siti N. Yaakub, University of Plymouth, 7 Sep 2022
arguments
    subj_id char
    t1_filename char
    ct_filename char
    output_dir char
    focus_coords (1,3) {mustBeNumeric}
    bowl_coords (1,3) {mustBeNumeric}
    focus_depth (1,1) {mustBeNumeric}
    transducer char
end
arguments (Repeating)
    varargin
end

% defaults
run_acoustic_sim = false;
run_thermal_sim = false;
pulse_length = 20e-3;	% pulse length [s]
pulse_rep_freq = 5;     % pulse repetition frequency [Hz]
stim_dur = 80;
run_cpp = 'matlab';

% specify transducer
[source_roc,diameters,freq] = get_transducer_specs(transducer);
[pressure,phase] = get_driving_params(focus_depth,transducer);

% replace with user-defined inputs, if given
if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'RunAcousticSim'
                run_acoustic_sim = varargin{arg_idx+1};
            case 'RunThermalSim'
                run_thermal_sim = varargin{arg_idx+1};
            case 'PulseLength'
                pulse_length = varargin{arg_idx+1};
            case 'PulseRepFreq'
                pulse_rep_freq = varargin{arg_idx+1};
            case 'StimDuration'
                stim_dur = varargin{arg_idx+1};
            case 'SourcePressure'
                pressure = varargin{arg_idx+1};
            case 'SourcePhase'
                phase = varargin{arg_idx+1};
            case 'RunCppCode'
                run_cpp = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

%% Simulation parameters
% Please change these if not using 500 kHz transducer.

% medium parameters
c_min               = 1500;     % sound speed [m/s]
c_max               = 3100;     % max. speed of sound in skull (F. A. Duck, 2013.) [m/s]
rho_min             = 1000;     % density [kg/m^3]
rho_max             = 1900;     % max. skull density [kg/m3]
alpha_power         = 1.43;     % Robertson et al., PMB 2017 usually between 1 and 3? from Treeby paper
alpha_coeff_water   = 0;        % [dB/(MHz^y cm)] 0.05 from Fomenko et al., 2020?
alpha_coeff_min     = 4;        %
alpha_coeff_max     = 8.7;      % [dB/(MHz cm)] Fry 1978 at 0.5MHz: 1 Np/cm (8.7 dB/cm) for both diploe and outer tables

% computational parameters
% ppw             = 3;        % number of points per wavelength
record_periods  = 3;        % number of periods to record
cfl             = 0.3;      % CFL number

%% Prepare skull & simulation, check model
%%% Skull preparation
hu_min 	= 300;	% minimum Hounsfield Unit in CT image
hu_max 	= 2000;	% maximum Hounsfield Unit in CT image

% Load CT image (nifti format)
% voxel size = 1 x 1 x 1 mm3, matrix size: varies, mostly 176x256x256
input_ct = niftiread(ct_filename);
t1_img = niftiread(t1_filename);
header = niftiinfo(ct_filename);
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
idx1 = zeros(3,2);
for ii = 1:3
    if shift_idx(ii,1)-100 < 0
        idx1(ii,1) = 1;
    elseif shift_idx(ii,1)-100 > 0
        idx1(ii,1) = shift_idx(ii,1)-100;
    end
    if shift_idx(ii,2)-100 <= size(input_ct,ii)
        idx1(ii,2) = shift_idx(ii,2)-100;
    elseif shift_idx(ii,2)-100 > size(input_ct,ii)
        idx1(ii,2) = size(input_ct,ii);
    end
end
idx2 = zeros(3,2);
idx2(:,1) = [shift_x, shift_y, shift_z] + [1,1,1];
idx2(:,2) = size(input_ct) + [shift_x, shift_y, shift_z];
for ii = 1:3
    if idx2(ii,1) < 0
        idx2(ii,1) = 1;
    end
    if idx2(ii,2) > size(model,ii)
        idx2(ii,2) = size(model,ii);
    end
end

new_t1(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2)) = ...
    t1_img(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2), idx1(3,1):idx1(3,2));
t1_img = new_t1;
clear new_t1;

%%% Medium properties
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

%%% Grid properties
% % calculate the grid spacing based on the PPW and F0
% dx = c_min / (ppw * freq);   % [m]
% get grid spacing from image header
dx = header.PixelDimensions(1)*1e-3; % [m]
% calculate ppw based on grid/image spacing and F0
ppw = c_min / (dx * freq);

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
ppp = round(ppw / cfl);
disp(['PPW = ' num2str(ppw)]);
disp(['CFL = ' num2str(cfl)]);

%%% Set source
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

if ~run_acoustic_sim
    % visualise transducer (check it is within the grid)
    model_mask = model;
    model_mask(model_mask>0)=1;
    tmp = source.p_mask+model_mask;
    volumeViewer(tmp);
    clear model_mask tmp;
    return
end

%% If 3D viewer of model looks OK, run simulation
% check if pressure and phase are empty
if isempty(pressure) || isempty(phase)
    error('No pressure or phase information provided.')
end
% driving parameters
source_amp = repmat(pressure, [1,4]);   % source pressure [Pa]
source_phase = deg2rad(phase);   % phase [rad]

% create time varying source
source_sig = createCWSignals(kgrid.t_array, freq, source_amp, source_phase);

% assign source signals (this step takes a long time)
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

%%% Set Sensor
% % set whole grid mask:
% sensor.mask = ones(Nx, Ny, Nz);
% set sensor mask within head mask: based on t1_img > 0
sensor.mask = zeros(size(model,1), size(model,2), size(model,3));
sensor.mask(t1_img>0) = 1;

% record the pressure
sensor.record = {'p'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * ppp + 1;

%%% Run simulation
% set input options
input_args = {'PMLSize', 10, 'PMLInside', true, 'PlotPML', true, ...
    'DisplayMask', source.p_mask};

switch run_cpp
    case 'matlab'
        % run C++ in Matlab
        sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
    case 'terminal'
        % run C++ in terminal
        % input and output filenames (these must have the .h5 extension)
        input_filename  = fullfile(output_dir, ['kwave_input_' subj_id '.h5']);
        output_filename = fullfile(output_dir, ['kwave_output_' subj_id '.h5']);
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', input_filename);
        % display the required syntax to run the C++ simulation
        disp(['Using a terminal window, navigate to the ' filesep 'binaries folder of the k-Wave Toolbox']);
        disp('Then, use the syntax shown below to run the simulation:');
        disp(['./kspaceFirstOrder-OMP -i ' input_filename ' -o ' output_filename ' --p_raw -s ' num2str(sensor.record_start_index)]);
        return
    case 'linux_system'
        % run C++ using system commands from matlab
        % define names of temporary input and output files
        input_filename = fullfile(output_dir, ['kwave_input_' subj_id '.h5']);
        output_filename = fullfile(output_dir, ['kwave_output_' subj_id '.h5']);
        
        % create input file
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', input_filename);
        
        % run simulation using bash commands
        kspacefct_path = fullfile(fileparts(which('kspaceFirstOrder3D')), ...
            'binaries', 'kspaceFirstOrder-OMP');
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

%%% Calculate pressure
% extract amplitude from the sensor data
p = extractAmpPhase(sensor_data.p, 1/kgrid.dt, freq, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
% % for whole grid mask:
% p = reshape(p, Nx, Ny, Nz);
% for sensor mask within head mask: based on t1_img > 0
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

p_focus = p(focus_coords(1), focus_coords(2), focus_coords(3));
isppa_focus = p_focus^2 / (2 * rho_min * c_min) * 1e-4;

%%% Check whether max pressure point is in brain
if norm(bowl_coords-[mx,my,mz])*dx*1e3 < focus_depth
    warning(['Maximum pressure point is not close to intended focus. ' ...
        'It is likely that the maximum pressure is found at the skull interface. ' ...
        'Attempting to adjust search to pressure recorded within the brain only... ' ...
        'Please check output!'])

    zz = mz-10;
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
end

%%% Summary
% Print summary to command window
disp(['Simulation output for: ' subj_id])
disp(['Shift = [' num2str(shift_x) ', ' num2str(shift_y) ', ' num2str(shift_z) ']'])
disp(['PPW = ' num2str(ppw)])
disp(['CFL = ' num2str(cfl)])
disp(['Coordinates of max pressure: [' num2str(mx-shift_x) ', ' num2str(my-shift_y) ', ' num2str(mz-shift_z) ']'])
disp(['Distance from transducer rear surface: ' num2str(norm(bowl_coords-[mx,my,mz])*dx*1e3) ' mm'])
disp(['Max Pressure = ' num2str(max_pressure * 1e-6) ' MPa'])
disp(['MI = ' num2str(MI)])
disp(['Isppa = ' num2str(Isppa) ' W/cm2'])
disp(['Ispta = ' num2str(Ispta*1e3) ' mW/cm2'])
disp(['Pressure at focus = ' num2str(p_focus * 1e-6) ' MPa']);
disp(['Isppa at focus = ' num2str(isppa_focus) ' W/cm2'])
disp(['-6dB focal volume = ' num2str(focal_vol) ' mm3'])
disp(' ')

% save output file in output_dir
clear source;
save(fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '.mat']));

% Write summary to spreadsheet
% create result file if it does not exist and write header
if ~exist(fullfile(output_dir, 'simulation_results.csv'), 'file')
 disp('Result file does not exist, creating file.')
 fileID = fopen(fullfile(output_dir, 'simulation_results.csv'), 'w' );
 fprintf(fileID, '%s\n', ...
    ['id, Focus coordinates, Transducer base coordinates, Focus depth, ' ...
    'PPW, CFL, Coordinates of max pressure, Distance from transducer rear surface (mm), ' ...
    'Max Pressure (MPa), MI, Isppa (W/cm2),' ...
    'Pressure at focus (MPa), Isppa at focus (W/cm2), ' ...
    '-6dB focal volume (mm3)']);
 fclose(fileID);
end

% write values
fileID = fopen(fullfile(output_dir, 'simulation_results.csv'),'a');
fprintf(fileID,['%s, %d %d %d, %d %d %d, %d, ' ...
    '%f, %f, %d %d %d, %f, ' ...
    '%f, %f, %f, %f, %f, %f\n'], ...
    subj_id, focus_coords, bowl_coords, focus_depth,  ...
    ppw, cfl, mx-shift_x, my-shift_y, mz-shift_z, norm(bowl_coords-[mx,my,mz])*dx*1e3, ...
    max_pressure * 1e-6, MI, Isppa, p_focus * 1e-6, isppa_focus, focal_vol);
fclose(fileID);

%%% Create Plots
figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(mx,:,:)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(p(mx,:,:))*1e-6,90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude')
saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_sag.jpg']));

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(mx,:,:)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(p(mx,:,:)>(0.5*max_pressure))*1e-6,90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'50% Acoustic Pressure Amplitude')
saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_sag_50%.jpg']));

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(:,my,:)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(p(:,my,:))*1e-6,90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude')
saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_cor.jpg']));

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(:,:,mz)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(p(:,:,mz))*1e-6,90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude')
saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_ax.jpg']));

%%% Save pressure field and -6dB volume as nifti file
p_out = zeros(size(input_ct),'double');
p_out(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2), idx1(3,1):idx1(3,2)) = ...
    p(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2));

header.Filename=[]; header.Filemoddate=[]; header.Filesize=[]; header.raw=[];
header.Datatype='double'; header.BitsPerPixel=32;
niftiwrite(p_out, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_pressurefield.nii']), header);

tmp_focusbin = p > 0.5*max_pressure;
tmp_focusbin = int16(tmp_focusbin);
focal_vol_bin = zeros(size(input_ct),'int16');
focal_vol_bin(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2), idx1(3,1):idx1(3,2)) = ...
    tmp_focusbin(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2));
clear tmp_focusbin;
header.Datatype='int16'; header.BitsPerPixel=16;
niftiwrite(focal_vol_bin, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_-6dBvolume.nii']), header);

%% Thermal Simulation
if ~run_thermal_sim
   return
end

% convert the absorption coefficient to nepers/m
alpha_np = db2neper(medium.alpha_coeff, medium.alpha_power) * ...
    (2 * pi * freq).^medium.alpha_power;

% reshape the data, and calculate the volume rate of heat deposition
Q = alpha_np .* p.^2 ./ (medium.density .* medium.sound_speed);

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

% set time step size
dt = on_time/1;

maxT1 = zeros(size(model));
for nb = 1:num_bursts
    disp(['burst #' num2str(nb)])
    % turn on heat source
    kdiff.Q = Q;
    % take time steps
    kdiff.takeTimeStep(round(on_time / dt), dt);

    % store the current temperature field
    T1 = kdiff.T;

    if max(T1(:)) > max(maxT1(:))
        maxT1 = T1;
    end

    % turn off heat source and take time steps
    kdiff.Q = 0;
    kdiff.takeTimeStep(round(off_time / dt), dt);

    % store the current temperature field
    T2 = kdiff.T;
end

[max_temp, idx_temp] = max(maxT1(:));
disp(['Max temperature rise = ' num2str(max_temp-source2.T0) ' °C'])
disp(['End max. temperature = ' num2str(max(T2(:))) ' °C'])
save(fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_thermalsim.mat']));

%%% Create thermal plots
[tx, ty, tz] = ind2sub(size(maxT1), idx_temp); % plot at max T coords

% plot the temperature maps at max temperature
figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(tx,:,:)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(maxT1(tx,:,:)),90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '°C');
title(ax1,'Max. temperature')
saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_therm_sag.jpg']));

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(:,ty,:)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(maxT1(:,ty,:)),90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '°C');
title(ax1,'Max. temperature')
saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_therm_cor.jpg']));

figure;
ax1 = axes;
imagesc(ax1, imrotate(squeeze(t1_img(:,:,tz)),90), [50,500]);
hold all;
ax2 = axes;
im2 = imagesc(ax2, imrotate(squeeze(maxT1(:,:,tz)),90));
im2.AlphaData = 0.5;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '°C');
title(ax1,'Max. temperature')
saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_therm_ax.jpg']));

% plot at max pressure
figure;
% plot the volume rate of heat deposition (yz)
subplot(3, 2, 1);
imagesc(kgrid.y_vec*1e3, kgrid.z_vec*1e3, imrotate(squeeze(Q(tx,:,:)*1e-7),90));
h = colorbar;
xlabel(h, '[kW/cm^2]');
ylabel('z-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Volume Rate Of Heat Deposition');

% plot the maximum temperature (yz)
subplot(3, 2, 2);
imagesc(kgrid.y_vec*1e3, kgrid.z_vec*1e3, imrotate(squeeze(maxT1(tx,:,:)),90));
h = colorbar;
xlabel(h, '[°C]');
ylabel('z-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Max. Temperature');

% plot the volume rate of heat deposition (xz)
subplot(3, 2, 3);
imagesc(kgrid.x_vec*1e3, kgrid.z_vec*1e3, imrotate(squeeze(Q(:,ty,:)*1e-7),90));
h = colorbar;
xlabel(h, '[kW/cm^2]');
ylabel('z-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Volume Rate Of Heat Deposition');

% plot the temperature after heating (xz)
subplot(3, 2, 4);
imagesc(kgrid.x_vec*1e3, kgrid.z_vec*1e3, imrotate(squeeze(maxT1(:,ty,:)),90));
h = colorbar;
xlabel(h, '[°C]');
ylabel('z-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Max. Temperature After Heating');

% plot the volume rate of heat deposition (xy)
subplot(3, 2, 5);
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, imrotate(Q(:,:,tz)*1e-7,90));
h = colorbar;
xlabel(h, '[kW/cm^2]');
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Volume Rate Of Heat Deposition');

% plot the temperature after heating (xy)
subplot(3, 2, 6);
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, imrotate(maxT1(:,:,tz),90));
h = colorbar;
xlabel(h, '[°C]');
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Max. Temperature After Heating');

% set colormap and enlarge figure window
colormap('turbo');
scaleFig(2, 3);
saveas(gcf, fullfile(output_dir, [subj_id '_tussim_skull_3D_' transducer '_thermalplots.jpg']));
