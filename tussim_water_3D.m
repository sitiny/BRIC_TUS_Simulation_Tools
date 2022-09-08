function tussim_water_3D(transducer, pressure, phase, varargin)
%TUSSIM_WATER_3D Run 3D k-wave acoustic simulation in free field for a 
%   given transducer, source amplitude and source phase.
% 
% Usage: 
%   tussim_water_3D(transducer, pressure, phase)
%   tussim_water_3D(transducer, pressure, phase, Name, Value)
%   tussim_water_3D('CTX500', 51590.0, [0, 319, 278, 237], ...
%       'RunThermalSim', true)
% 
% Inputs:
%   transducer: Options are 'CTX500', 'CTX560' or 'CTX250'.
%   pressure:   Source pressure in Pa.
%   phase:      4-element array of phases of each transducer element in 
%               degrees for the focal depth required.
%                 
% Optional Name-Value pair inputs:
%   'RunThermalSim':    Boolean controlling whether to run thermal 
%                       simulation(default: false).
%   'PulseLength':      Pulse length in s (default: 20e-3 s).
%   'PulseRepFreq':     Pulse repetition frequency in Hz (default: 5 Hz).
%   'StimDuration':     Duration of TUS in s (default: 80 s).
% 
% WARNING: Default acoustic values in function are for 500 kHz transducers. 
% Please set your own values if using central frequency other than 500kHz.
% 
% Dependencies:
%   k-Wave Toolbox (http://www.k-wave.org)
%   kArray (http://www.k-wave.org/downloads/kWaveArray_alpha_0.3.zip)
%     	Copyright (C) 2009-2017 Bradley Treeby
% 
% Author: Siti N. Yaakub, University of Plymouth, 7 Sep 2022
arguments
    transducer char
    pressure (1,1) {mustBeNumeric}
    phase (1,4) {mustBeNumeric}
end
arguments (Repeating)
    varargin
end

% defaults
run_thermal_sim = false;
pulse_length = 20e-3;	% pulse length [s]
pulse_rep_freq = 5;     % pulse repetition frequency [Hz]
stim_dur = 80;

% replace with user-defined inputs, if given
if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'RunThermalSim'
                run_thermal_sim = varargin{arg_idx+1};
            case 'PulseLength'
                pulse_length = varargin{arg_idx+1};
            case 'PulseRepFreq'
                pulse_rep_freq = varargin{arg_idx+1};
            case 'StimDuration'
                stim_dur = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

% transducer specs
[source_roc,diameters,freq] = get_transducer_specs(transducer);

% driving parameters
source_amp = repmat(pressure, [1,4]); % source pressure [Pa]
source_phase = deg2rad(phase); % phase [rad]

% source parameters
freq = freq*1e3; % source frequency [Hz]
bx = 10; fx = 54;

% medium parameters
c0              = 1500;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]
alpha_power     = 1.43;     % Robertson et al., PMB 2017
alpha_coeff     = 0.05;     % [dB/(MHz^y cm)] Fomenko et al., 2020

% grid parameters
axial_size      = 128e-3;    % total grid size in the axial dimension [m]
lateral_size    = 128e-3;    % total grid size in the lateral dimension [m]

% computational parameters
ppw             = 3;        % number of points per wavelength
record_periods  = 3;        % number of periods to record
cfl             = 0.3;      % CFL number

%% prepare simulation
% calculate the grid spacing based on the PPW and F0
dx = c0 / (ppw * freq);   % [m]

% compute the size of the grid
Nx = roundEven(axial_size / dx);
Ny = roundEven(lateral_size / dx);
Nz = Ny;

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% compute points per temporal period
ppp = round(ppw / cfl);

% compute corresponding time spacing
dt = 1 / (ppp * freq);

% calculate the number of time steps to reach steady state
t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2) / c0; 

% create the time array using an integer number of points per period
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
disp(['PPW = ' num2str(c0 / (dx * freq))]);
disp(['CFL = ' num2str(c0 * dt / dx)]);

% create time varying source
source_sig = createCWSignals(kgrid.t_array, freq, source_amp, source_phase);

% set bowl position and orientation
bowl_pos = [kgrid.x_vec(bx), 0, 0];
focus_pos = [kgrid.x_vec(fx), 0, 0];

% create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.01, 'UpsamplingRate', 16);

% add bowl shaped element
karray.addAnnularArray(bowl_pos, source_roc, diameters, focus_pos)

% assign binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% assign source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

% assign medium properties
medium.sound_speed = c0;
medium.density = rho0;
medium.alpha_power = alpha_power;
medium.alpha_coeff = alpha_coeff;

% set sensor mask 
sensor.mask = ones(Nx, Ny, Nz);

% record the pressure
sensor.record = {'p'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * ppp + 1;

%% Run acoustic simulation
% set input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', true, 'DisplayMask', 'off'};
sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});

% extract amplitude from the sensor data
p = extractAmpPhase(sensor_data.p, 1/kgrid.dt, freq, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
p = reshape(p, Nx, Ny, Nz);

% extract pressure on axis
amp_on_axis = p(:, Ny/2 + 1, Nz/2 + 1);

% define axis vectors for plotting
x_vec = kgrid.x_vec - kgrid.x_vec(1);
y_vec = kgrid.y_vec;

figure;
plot(x_vec*1e3, amp_on_axis*1e-6, 'b.');
set(gca, 'XLim', [0, Nx*dx] * 1e3);
xlabel('Axial Position [mm]');
ylabel('Pressure [MPa]');
title('Axial Pressure');

% plot the pressure field 
figure;
imagesc(y_vec*1e3, x_vec*1e3, imrotate(p(:,:,Nz/2+1),90));
colormap('turbo');
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
axis image;
title('Pressure Field');

% calculate acoustic intensities
[max_pressure, idx] = max(p(:)); % [Pa]
[mx, my, mz] = ind2sub(size(p), idx);

Isppa = max_pressure^2 / (2 * medium.density * medium.sound_speed); % [W/m2]
Isppa = Isppa * 1e-4; % [W/cm2]
Ispta = Isppa * pulse_length * pulse_rep_freq; % [W/cm2]

% MI = max_pressure (in MPa) / sqrt freq (in MHz)
MI = max_pressure * 1e-6 / sqrt(freq * 1e-6);

disp(['Max Pressure = ' num2str(max_pressure * 1e-6) ' MPa'])
disp(['MI = ' num2str(MI)])
disp(['Coordinates of max pressure: ' num2str(mx) ', ' num2str(my) ', ' num2str(mz)])
disp(['Distance from transducer rear surface: ' num2str(norm([bx,find(kgrid.x_vec==0),find(kgrid.x_vec==0)]-[mx,my,mz])*dx*1e3) ' mm'])
disp(['Isppa = ' num2str(Isppa) ' W/cm2'])
disp(['Ispta = ' num2str(Ispta*1e3) ' mW/cm2'])

%% Run thermal simulation
if ~run_thermal_sim
   return 
end

% convert the absorption coefficient to nepers/m
alpha_np = db2neper(medium.alpha_coeff, medium.alpha_power) * ...
    (2 * pi * freq).^medium.alpha_power;

% reshape the data, and calculate the volume rate of heat deposition
Q = alpha_np .* p.^2 ./ (medium.density .* medium.sound_speed);

% clear the input structures
clear medium source sensor;

% set the background temperature and heating term
source.Q = Q;
source.T0 = 37;

% define medium properties related to diffusion
medium.density              = 1020;     % [kg/m^3]
medium.thermal_conductivity = 0.5;      % [W/(m.K)]
medium.specific_heat        = 3600;     % [J/(kg.K)]

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium, source, [], 'DisplayUpdates', false, 'PlotSim', false);

% set source on time and off time
on_time  = pulse_length;  % [s]
off_time = 1/pulse_rep_freq - pulse_length;  % [s]
num_bursts = stim_dur * pulse_rep_freq; % 

% set time step size
dt = on_time/1;

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

% plot the thermal dose and lesion map
figure;
% plot the acoustic pressure
subplot(2, 2, 1);
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, imrotate(p(:,:,Nz/2+1) * 1e-6,90));
h = colorbar;
xlabel(h, '[MPa]');
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Acoustic Pressure Amplitude');

% plot the volume rate of heat deposition
subplot(2, 2, 2);
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, imrotate(Q(:,:,Nz/2+1) * 1e-7,90));
h = colorbar;
xlabel(h, '[kW/cm^2]');
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Volume Rate Of Heat Deposition');

% plot the temperature after heating
subplot(2, 2, 3);
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, imrotate(T1(:,:,Nz/2+1),90));
h = colorbar;
xlabel(h, '[degC]');
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
title('Temperature After Heating');

% plot the temperature after cooling
subplot(2, 2, 4);
imagesc(kgrid.x_vec * 1e3, kgrid.y_vec * 1e3, T2(:,:,Nz/2+1).');%, [37, 37.1]);
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

% --------------------
% SUMMARY
% --------------------

disp(['Max Pressure = ' num2str(max_pressure * 1e-6) ' MPa'])
disp(['MI = ' num2str(MI)])
disp(['Coordinates of max pressure: ' num2str(mx) ', ' num2str(my) ', ' num2str(mz)])
disp(['Distance from transducer rear surface: ' num2str(norm([bx,find(kgrid.x_vec==0),find(kgrid.x_vec==0)]-[mx,my,mz])*dx*1e3) ' mm'])
disp(['Isppa = ' num2str(Isppa) ' W/cm2'])
disp(['Ispta = ' num2str(Ispta*1e3) ' mW/cm2'])
disp(['Max temperature rise = ' num2str(max(T1(:))) ' degC'])