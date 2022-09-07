% parameters to change depending on the transducer and pulse protocol
% pulse parameters
pulse_length    = 20e-3;	% pulse length [s]
pulse_rep_freq  = 5;        % pulse repetition frequency [Hz]
stim_dur        = 80;       % duration of stimulation [s]

% source parameters
freq        = 500e3;            % source frequency [Hz]
source_roc  = 63.2e-3;          % bowl radius of curvature [m]

% aperture diameters of the elements given an inner, outer pairs [m]
diameters       = [0 1.28; 1.3 1.802; 1.82 2.19; 2.208 2.52].' * 0.0254;

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
ppw             = 3;        % number of points per wavelength
record_periods  = 3;        % number of periods to record
cfl             = 0.3;      % CFL number

save simulation_params.mat