function tussim_skull_3D(t1_filename, ct_filename, output_dir, ...
    focus_coords_in, bowl_coords_in, focus_depth, varargin)
%TUSSIM_SKULL_3D Run 3D k-wave acoustic simulation transcranially with skull
%   estimated using CT (or pseudo-CT) images.
%
% The transducer is modelled on the NeuroFUS CTX-500-4 with free-field Isppa
% of 20 W/cm^2. If you would like the simulation for a different free-field
% Isppa or a different 4-element transducer, you will have to provide your
% own source pressure and phase as optional Name-Value paired input arguments.
%
% The simulation grid is 256x256x256, with grid spacing determined by PPW
% (at recommended PPW = 6, grid spacing = 0.5x0.5x0.5mm^3).
%
% You will need to supply a co-registered T1-weighted MR image and CT (or
% pseudo-CT) image for use in the simulations. These images are recommended
% to be resampled to 1mm^3 isotropic voxels.
%
% Running the script without the acoustic or thermal simulation allows you
% to check that the transducer position relative to the head is correct.
% If you are satisfied, re-run the script with ('RunAcousticSim', true).
%
% Usage:
%   tussim_skull_3D(t1_filename, ct_filename, output_dir, ...
%       focus_coords, bowl_coords, focus_depth)
%   tussim_skull_3D(t1_filename, ct_filename, output_dir, ...
%       focus_coords, bowl_coords, focus_depth, Name, Value)
%   tussim_skull_3D('sub-test01_t1w.nii', 'sub-test01_pct.nii', 'output_dir', ...
%       [99, 161, 202], [90, 193, 262], 60, 'RunAcousticSim', true);
%
% Inputs:
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
%
% Optional Name-Value pair inputs:
%   'PPW':              Points per wavelength (default: 3, recommended: 6).
%   'RunAcousticSim':   Boolean controlling whether to run acoustic
%                       simulation (default: false).
%   'RunThermalSim':    Boolean controlling whether to run thermal
%                       simulation(default: false).
%   'PulseDur':         Pulse duration in s (default: 0.02 s).
%   'PulseRepInt':      Pulse repetition interval in s (default: 0.2s).
%   'PulseTrainDur':    Pulse train duration in s (default: 80 s).
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
%   k-Wave Toolbox v 1.4 (http://www.k-wave.org)
%     	Copyright (C) 2009-2017 Bradley Treeby
%
% Author: Siti N. Yaakub, University of Plymouth, 7 Sep 2022
%         (edited 2 Jun 2023)
arguments
    t1_filename char
    ct_filename char
    output_dir char
    focus_coords_in (1,3) {mustBeNumeric}
    bowl_coords_in (1,3) {mustBeNumeric}
    focus_depth (1,1) {mustBeNumeric}
end
arguments (Repeating)
    varargin
end

% defaults
run_acoustic_sim = false;
run_thermal_sim = false;
run_cpp = 'linux_system';
ppw = 3;

% replace with user-defined inputs, if given
if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'PPW'
                ppw = varargin{arg_idx+1};
            case 'RunAcousticSim'
                run_acoustic_sim = varargin{arg_idx+1};
            case 'RunThermalSim'
                run_thermal_sim = varargin{arg_idx+1};
            case 'PulseDur'
                pulse_dur = varargin{arg_idx+1};
            case 'PulseRepInt'
                pulse_rep_int = varargin{arg_idx+1};
            case 'PulseTrainDur'
                pulse_train_dur = varargin{arg_idx+1};
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

% specify transducer
[source_roc,diameters,freq] = get_transducer_specs('CTX500');
if ~exist('pressure', 'var') || ~exist('phase', 'var')
    disp('Source pressure and phases not set, using default values for focus depth.')
    [pressure,phase] = get_driving_params(focus_depth,'CTX500');
end

%% Simulation parameters
% Please change these if not using 500 kHz transducer.

% medium parameters
c_min               = 1500;     % sound speed [m/s]
c_max               = 3100;     % max. speed of sound in skull (F. A. Duck, 2013.) [m/s]
rho_min             = 1000;     % density [kg/m^3]
rho_max             = 1900;     % max. skull density [kg/m3]
alpha_power         = 1.43;     % Robertson et al., PMB 2017 usually between 1 and 3? from Treeby paper
alpha_coeff_water   = 0;        % [dB/(MHz^y cm)] close to 0 (Mueller et al., 2017), see also 0.05 Fomenko et al., 2020?
alpha_coeff_min     = 4;        %
alpha_coeff_max     = 8.7;      % [dB/(MHz cm)] Fry 1978 at 0.5MHz: 1 Np/cm (8.7 dB/cm) for both diploe and outer tables

% computational parameters
record_periods  = 3;        % number of periods to record
cfl             = 0.3;      % CFL number

%% Prepare skull & simulation, check model
%%% Skull preparation
hu_min 	= 300;	% minimum skull HU
hu_max 	= 2000;	% maximum skull HU

% Load CT image (nifti format)
% voxel size = 1 x 1 x 1 mm3, matrix size: varies
input_ct = niftiread(ct_filename);
input_t1 = niftiread(t1_filename);
header = niftiinfo(ct_filename);

% calculate the grid spacing based on PPW and F0
dx = c_min / (ppw * freq); % in mm

% resample input images to grid res (iso)
scale_factor = round(header.PixelDimensions/(dx*1e3),2);
ct_img = imresize3(input_ct, 'cubic', 'Scale', scale_factor);
t1_img = imresize3(input_t1, 'cubic', 'Scale', scale_factor);

focus_coords = round(focus_coords_in.*scale_factor);
bowl_coords = round(bowl_coords_in.*scale_factor);

% update hu_max
ct_max = max(ct_img(:));
if ct_max < hu_max
    hu_max = ct_max;
end
clear ct_max;

% truncate CT HU (see Marsac et al., 2017)
skull_model = ct_img;
skull_model(skull_model < hu_min) = 0; % only use HU for skull acoustic properties
skull_model(skull_model > hu_max) = hu_max;

% centre grid at midpoint between transducer and focus coordinates
midpoint = round((bowl_coords + focus_coords)/2);

% pad images by 128 on each side
padx = 128;
tmp_model = zeros(size(skull_model,1)+padx*2, size(skull_model,2)+padx*2, size(skull_model,3)+padx*2);
tmp_model(padx+1:size(skull_model,1)+padx, ...
    padx+1:size(skull_model,2)+padx, ...
    padx+1:size(skull_model,3)+padx) = skull_model;
tmp_midpoint = midpoint+padx;

% centre on midpoint
% grid size = 256x256x256, new midpoint coords = [128,128,128]
shift_idx = [tmp_midpoint(1)-padx+1,tmp_midpoint(1)+padx;
    tmp_midpoint(2)-padx+1,tmp_midpoint(2)+padx; ...
    tmp_midpoint(3)-padx+1,tmp_midpoint(3)+padx];
model = tmp_model(shift_idx(1,1):shift_idx(1,2), ...
    shift_idx(2,1):shift_idx(2,2), ...
    shift_idx(3,1):shift_idx(3,2));

shift_x = padx-midpoint(1);
shift_y = padx-midpoint(2);
shift_z = padx-midpoint(3);

bowl_coords     = bowl_coords + [shift_x, shift_y, shift_z];	% centre of rear surface of transducer [grid points]
focus_coords    = focus_coords + [shift_x, shift_y, shift_z];  % point on the beam axis of the transducer [grid points]

% move t1 image
new_t1 = zeros(size(model));
idx1 = zeros(3,2);
for ii = 1:3
    if shift_idx(ii,1)-padx <= 0
        idx1(ii,1) = 1;
    elseif shift_idx(ii,1)-padx > 0
        idx1(ii,1) = shift_idx(ii,1)-padx;
    end
    if shift_idx(ii,2)-padx <= size(ct_img,ii)
        idx1(ii,2) = shift_idx(ii,2)-padx;
    elseif shift_idx(ii,2)-padx > size(ct_img,ii)
        idx1(ii,2) = size(ct_img,ii);
    end
end
idx2 = zeros(3,2);
idx2(:,1) = [shift_x, shift_y, shift_z] + [1,1,1];
idx2(:,2) = size(ct_img) + [shift_x, shift_y, shift_z];
for ii = 1:3
    if idx2(ii,1) <= 0
        idx2(ii,1) = 1;
    end
    if idx2(ii,2) > size(model,ii)
        idx2(ii,2) = size(model,ii);
    end
end

new_t1(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2)) = ...
    t1_img(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2), idx1(3,1):idx1(3,2));
t1_img = new_t1;
clear new_t1 tmp_model;

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

% calculate the actual CFL after adjusting for dt
cfl = c_min * dt / dx;
ppp = round(ppw / cfl);
disp(['PPW = ' num2str(ppw)]);
disp(['CFL = ' num2str(cfl)]);
disp(['PPP = ' num2str(ppp)]);

%%% Set source
% set bowl position and orientation
bowl_pos = [kgrid.x_vec(bowl_coords(1)), ...
    kgrid.y_vec(bowl_coords(2)), ...
    kgrid.z_vec(bowl_coords(3))];
focus_pos = [kgrid.x_vec(focus_coords(1)), ...
    kgrid.y_vec(focus_coords(2)), ...
    kgrid.z_vec(focus_coords(3))];

% create empty kWaveArray
karray = kWaveArray('SinglePrecision', true);

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
if run_acoustic_sim
    
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
    sensor.mask = ones(Nx, Ny, Nz);
    %     % set sensor mask within head mask: based on t1_img > 0
    %     sensor.mask(t1_img>0) = 1;
    
    % record the pressure
    sensor.record = {'p'};
    
    % average only the final few periods when the field is in steady state
    sensor.record_start_index = kgrid.Nt - record_periods * ppp + 1;
    
    %%% Run simulation
    % set input options
    input_args = {'PMLSize', 10, 'PMLInside', true, 'PlotPML', true, ...
        'DisplayMask', source.p_mask};
    
    if ~exist(fullfile(output_dir), 'dir')
        mkdir(output_dir)
    end
    switch run_cpp
        case 'matlab'
            % run C++ in Matlab
            sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
        case 'terminal'
            % run C++ in terminal
            % input and output filenames (these must have the .h5 extension)
            input_filename  = fullfile(output_dir, 'kwave_input.h5');
            output_filename = fullfile(output_dir, 'kwave_output.h5');
            sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', input_filename);
            % display the required syntax to run the C++ simulation
            disp(['Using a terminal window, navigate to the ' filesep 'binaries folder of the k-Wave Toolbox']);
            disp('Then, use the syntax shown below to run the simulation:');
            disp(['./kspaceFirstOrder-OMP -i ' input_filename ' -o ' output_filename ' --p_raw -s ' num2str(sensor.record_start_index)]);
            return
        case 'linux_system'
            % run C++ using system commands from matlab
            % define names of temporary input and output files
            input_filename = fullfile(output_dir, 'kwave_input.h5');
            output_filename = fullfile(output_dir, 'kwave_output.h5');
            
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
            system(['rm ' input_filename ' ' output_filename]);
    end
    
    %%% Calculate pressure
    % extract amplitude from the sensor data
    p = extractAmpPhase(sensor_data.p, 1/kgrid.dt, freq, ...
        'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
    
    % reshape data
    % for whole grid mask:
    p = reshape(p, Nx, Ny, Nz);
    %     % for sensor mask within head mask: based on t1_img > 0
    %     tmp_p = p;
    %     p = zeros(size(model));
    %     p(t1_img>0) = tmp_p;
    %     clear tmp_p;
    
    % calculate acoustic intensities
    [max_pressure, idx] = max(p(:)); % [Pa]
    [mx, my, mz] = ind2sub(size(p), idx);
    
    if model(mx, my, mz) > 1
        Isppa = max_pressure^2 / (2 * max(medium.density(:)) * max(medium.sound_speed(:))); % [W/m2]
    elseif model(mx, my, mz) == 0
        Isppa = max_pressure^2 / (2 * rho_min * c_min); % [W/m2]
    end
    Isppa = Isppa * 1e-4; % [W/cm2]
    %     Ispta = Isppa * pulse_length * pulse_rep_freq; % [W/cm2]
    
    % MI = max_pressure (in MPa) / sqrt freq (in MHz)
    MI = max_pressure * 1e-6 / sqrt(freq * 1e-6);
    
    % find -6dB focal volume
    % get largest connected component - probably the main focus
    tmp_focal_vol = int16(p>0.5*max(p(:)));
    cc = bwconncomp(tmp_focal_vol);
    focal_vol = length(cc.PixelIdxList{1})*(dx*1e3)^3;
    clear tmp_focal_vol cc;
    
    p_focus = p(focus_coords(1), focus_coords(2), focus_coords(3));
    isppa_focus = p_focus^2 / (2 * rho_min * c_min) * 1e-4;
    
    %%% Check whether max pressure point is in brain
    if norm(focus_coords-[mx,my,mz])*dx*1e3 > 5
        warning(['Maximum pressure point is more than 5 mm away from the intended focus. ' ...
            'It is likely that the maximum pressure is at the skull interface. ' ...
            'Please check output!'])
    end
    
    %%% Create Plots
    figure;
    ax1 = axes; imagesc(ax1, imrotate(squeeze(t1_img(mx,:,:)),90), [50,500]);
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
    saveas(gcf, fullfile(output_dir, 'pressure_sag.jpg'));
    
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
    saveas(gcf, fullfile(output_dir, 'pressure_sag_50%.jpg'));
    
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
    saveas(gcf, fullfile(output_dir, 'pressure_cor.jpg'));
    
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
    saveas(gcf, fullfile(output_dir, 'pressure_ax.jpg'));
    
    %%% Save pressure field and -6dB volume as nifti file
    p_out = zeros(size(ct_img),'double');
    p_out(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2), idx1(3,1):idx1(3,2)) = ...
        p(idx2(1,1):idx2(1,2), idx2(2,1):idx2(2,2), idx2(3,1):idx2(3,2));
    p_out = imresize3(p_out, 'cubic', 'Scale', 1./scale_factor);
    header.ImageSize = size(p_out);
    header.Filename=[]; header.Filemoddate=[]; header.Filesize=[]; header.raw=[];
    header.Datatype='double'; header.BitsPerPixel=32;
    niftiwrite(p_out, fullfile(output_dir, 'pressure_field.nii'), header);
    
    focal_vol_bin = int16(p_out > 0.5*max_pressure);
    cc = bwconncomp(focal_vol_bin);
    focal_vol_lcc = int16(zeros(size(p_out)));
    focal_vol_lcc(cc.PixelIdxList{1}) = 1;
    header.Datatype='int16'; header.BitsPerPixel=16;
    niftiwrite(focal_vol_lcc, fullfile(output_dir, 'focal_volume_bin.nii'), header);
    
    % find max pressure point on original image
    [max_pressure, ~] = max(p_out(logical(focal_vol_lcc))); % [Pa]
    idx = find(p_out==max_pressure);
    [mx, my, mz] = ind2sub(size(p_out), idx);
    
    %%% Summary
    % Print summary to command window
    disp(['PPW = ' num2str(ppw)])
    disp(['CFL = ' num2str(cfl)])
    disp(['Coordinates of max pressure: [' num2str(mx) ', ' num2str(my) ', ' num2str(mz) ']'])
    disp(['Max Pressure = ' num2str(max_pressure * 1e-6) ' MPa'])
    disp(['MI = ' num2str(MI)])
    disp(['Isppa = ' num2str(Isppa) ' W/cm2'])
    disp(['Pressure at focus = ' num2str(p_focus * 1e-6) ' MPa']);
    disp(['Isppa at focus = ' num2str(isppa_focus) ' W/cm2'])
    disp(['-6dB focal volume = ' num2str(focal_vol) ' mm3'])
    disp(' ')
    
    % save output file in output_dir
    clear source sensor_data;
    save(fullfile(output_dir, 'acoustic_sim_output.mat'));
    
    % Write summary to spreadsheet
    % create result file if it does not exist and write header
    if ~exist(fullfile(output_dir, 'simulation_results.csv'), 'file')
        disp('Result file does not exist, creating file.')
        fileID = fopen(fullfile(output_dir, 'simulation_results.csv'), 'w' );
        fprintf(fileID, '%s\n', ...
            ['Output directory, Focus coordinates, Bowl coordinates, Focus depth, ' ...
            'PPW, CFL, PPP, Coordinates of max pressure, ' ...
            'Max Pressure (MPa), MI, Isppa (W/cm2),' ...
            'Pressure at focus (MPa), Isppa at focus (W/cm2), ' ...
            '-6dB focal volume (mm3)']);
        fclose(fileID);
    end
    
    % write values
    fileID = fopen(fullfile(output_dir, 'simulation_results.csv'),'a');
    fprintf(fileID,['%s, %d %d %d, %d %d %d, %d, ' ...
        '%f, %f, %f, %d %d %d, ' ...
        '%f, %f, %f, %f, %f, %f\n'], ...
        output_dir, focus_coords_in, bowl_coords_in, focus_depth,  ...
        ppw, cfl, ppp, mx, my, mz, ...
        max_pressure * 1e-6, MI, Isppa, p_focus * 1e-6, isppa_focus, focal_vol);
    fclose(fileID);
end
%% Thermal Simulation
if run_thermal_sim
    
    % check input pulse parameters
    if ~exist('pulse_dur', 'var')
        disp('Thermal simulation will be run with default pulse length.');
        pulse_dur = 20e-3;	% pulse length [s]
    end
    if ~exist('pulse_rep_int', 'var')
        disp('Thermal simulation will be run with default pulse repetition frequency.')
        pulse_rep_int = 0.2;     % pulse repetition interval [s]
    end
    if ~exist('pulse_train_dur', 'var')
        disp('Thermal simulation will be run with default pulse train duration.')
        pulse_train_dur = 80;
    end
    
    % load acoustic simulation
    load(fullfile(output_dir, 'acoustic_sim_output.mat'), 'medium', 'p');
    
    % convert the absorption coefficient to nepers/m
    alpha_np = db2neper(medium.alpha_coeff, medium.alpha_power) * ...
        (2 * pi * freq).^medium.alpha_power;
    
    % reshape the data, and calculate the volume rate of heat deposition
    Q = alpha_np .* p.^2 ./ (medium.density .* medium.sound_speed);
    
    % downsample everything so thermal simulation runs faster
    NTherm = round([Nx,Ny,Nz]./scale_factor);
    dxT= dx*scale_factor;
    kgridTherm = kWaveGrid(NTherm(1),dxT(1),NTherm(2),dxT(2),NTherm(3),dxT(3));
    
    Q = imresize3(Q, 'cubic', 'Scale', 1./scale_factor);
    
    % set the background temperature and heating term
    sourceTherm.Q = Q;
    sourceTherm.T0 = 37;
    
    % set sensor mask within skull mask + some dilation (assuming max T
    % rise will be at skull/brain interface)
    modelTherm = imresize3(model, 'cubic', 'Scale', 1./scale_factor);
    sensorTherm.mask = modelTherm;
    sensorTherm.mask(sensorTherm.mask>0) = 1;
    sensorTherm.mask = imdilate(sensorTherm.mask,strel('cube',20));
    
    % define medium properties related to diffusion
    % ref: https://itis.swiss/virtual-population/tissue-properties/database
    % in skull
    mediumTherm.density = rho_min + (rho_max - rho_min) * (modelTherm - 0) / (hu_max - 0);
    mediumTherm.thermal_conductivity = zeros(size(modelTherm));
    mediumTherm.thermal_conductivity(modelTherm > 0) = 0.32;
    mediumTherm.specific_heat = 826 + (2524 - 826) * (1-((modelTherm - hu_min) / (hu_max - hu_min)));
    % in water
    mediumTherm.density(modelTherm == 0)              = rho_min;     % [kg/m^3]
    mediumTherm.thermal_conductivity(modelTherm == 0) = 0.6;      % [W/(m.K)]
    mediumTherm.specific_heat(modelTherm == 0)        = 4178;     % [J/(kg.K)]
    
    % create kWaveDiffusion object
    kdiff = kWaveDiffusion(kgridTherm, mediumTherm, sourceTherm, sensorTherm.mask, 'DisplayUpdates', false, 'PlotSim', false);
    
    % set source on time and off time
    on_time  = pulse_dur;  % [s]
    off_time = pulse_rep_int - pulse_dur;  % [s]
    num_pulses = pulse_train_dur / pulse_rep_int;
    
    % set time step size
    dt = on_time/10;
    
    %     step_count = 0;
    maxT1 = zeros(size(model));
    idxT = 1;
    timepoint(idxT) = 0;
    stepT(idxT) = sourceTherm.T0;
    tic
    for np = 1:num_pulses
        disp(['pulse #' num2str(np)])
        
        % turn on heat source
        kdiff.Q = Q;
        % take time steps
        kdiff.takeTimeStep(round(on_time / dt), dt);
        % store the current temperature field
        T1 = kdiff.T;
        if max(T1(:)) > max(maxT1(:))
            maxT1 = T1;
        end
        idxT = idxT + 1;
        timepoint(idxT) = timepoint(idxT-1)+on_time;
        stepT(idxT) = max(T1(:));
        
        % turn off heat source and take time steps
        kdiff.Q = 0;
        kdiff.takeTimeStep(round(off_time / dt), dt);
        % store the current temperature field
        T2 = kdiff.T;
        idxT = idxT + 1;
        timepoint(idxT) = timepoint(idxT-1)+off_time;
        stepT(idxT) = max(T2(:));
    end
    toc
    [max_temp, idx_temp] = max(maxT1(:));
    disp(['Max temperature rise = ' num2str(max_temp-sourceTherm.T0) ' °C'])
    disp(['End max. temperature = ' num2str(max(T2(:))) ' °C'])
    save(fullfile(output_dir, 'thermal_sim_output.mat'));
    
    %%% Create thermal plots
    [tx, ty, tz] = ind2sub(size(maxT1), idx_temp); % plot at max T coords
    
    t1_imgT = imresize3(t1_img, 'cubic', 'Scale', 1./scale_factor);
    
    % plot the temperature maps at max temperature
    figure;
    ax1 = axes;
    imagesc(ax1, imrotate(squeeze(t1_imgT(tx,:,:)),90), [50,500]);
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
    saveas(gcf, fullfile(output_dir, 'temperature_sag.jpg'));
    
    figure;
    ax1 = axes;
    imagesc(ax1, imrotate(squeeze(t1_imgT(:,ty,:)),90), [50,500]);
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
    saveas(gcf, fullfile(output_dir, 'temperature_cor.jpg'));
    
    figure;
    ax1 = axes;
    imagesc(ax1, imrotate(squeeze(t1_imgT(:,:,tz)),90), [50,500]);
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
    saveas(gcf, fullfile(output_dir, 'temperature_ax.jpg'));
    
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
    saveas(gcf, fullfile(output_dir, 'thermal_plots.jpg'));
end
