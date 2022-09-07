% driving parameters (depends on focus depth)
idx = find(driving_params.dist == focus_depth);
source_amp      = repmat(driving_params.amp{idx}, [1,4]);   % source pressure [Pa]
source_phase    = deg2rad([driving_params.phase{idx,:}]);   % phase [rad]

% create time varying source
source_sig = createCWSignals(kgrid.t_array, freq, source_amp, source_phase);

% assign source signals (this step takes a long time)
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

% --------------------
% SENSOR
% --------------------

% set sensor mask (restrict mask to head mask based on t1_img - assumes t1_img has been masked for pct creation)
% sensor.mask = ones(Nx, Ny, Nz);
sensor.mask = zeros(size(model,1), size(model,2), size(model,3));
sensor.mask(t1_img>0) = 1;

% record the pressure
sensor.record = {'p'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * round(ppw/cfl) + 1;

%%% Run simulation

% --------------------
% ACOUSTIC SIMULATION
% --------------------

% set input options
input_args = {'PMLSize', 10, 'PMLInside', true, 'PlotPML', true, 'DisplayMask', source.p_mask};

% run C++ in Matlab
sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
toc;

% --------------------
% CALCULATE PRESSURE
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