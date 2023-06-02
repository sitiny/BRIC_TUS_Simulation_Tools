addpath('/Users/sitiyaakub/Documents/MATLAB/toolboxes/k-wave-toolbox-version-1.4/k-Wave')
addpath(genpath('/Users/sitiyaakub/Documents/Analysis/BRIC_TUS_Simulations/github/'))

addpath('~/Software/k-wave-toolbox-version-1.4/k-Wave') % copied binaries from v1.3
addpath(genpath('/usr/local/Software/BRIC_TUS_Simulation_Tools-main/utils'))

% subj_id         = 'IALB0252';      % subject ID
focus_coords_in = [98, 57, 120] +1;
bowl_coords_in = [123, 5, 118] +1;
focus_depth     = 46;                   % to nearest mm
% transducer = 'CTX500';

% base_dir = '/2TBdisk/Siti-Data/scotoma/sub-IALB0252';
base_dir = '/Users/sitiyaakub/Documents/Analysis/Scotoma/test_sims/sub-IALB0252';
output_dir = fullfile(base_dir, 'v1_ppw6_newfunc');
ct_filename = fullfile(base_dir, 'sub-IALB0252_T1w_n4_masked_pct.nii');
t1_filename = fullfile(base_dir, 'sub-IALB0252_T1w_n4_masked.nii');

% sim parameters
pulse_dur = 100e-3; 
pulse_rep_int = 1/5;
pulse_train_dur = 600e-3;
run_cpp = 'linux_system';
% isppa = 54;
pressure = 46550; %69017.41 for 54.49 Isppa
phase = [0 15.3 30.6 46]; % 46 mm

tussim_skull_3D(t1_filename, ct_filename, output_dir, ...
    focus_coords_in, bowl_coords_in, focus_depth, 'PPW', 6, ...
    'RunAcousticSim', true, 'RunThermalSim', true, 'PulseDur', pulse_dur, ...
    'PulseRepInt', pulse_rep_int, 'PulseTrainDur', pulse_train_dur, ...
    'SourcePressure', pressure, 'SourcePhase', phase)

%%
size0=size(t1_img);
figure;
slice(double(t1_img),size0(2)/2,size0(1)/2,size0(3)/2);
shading interp, colormap gray;
title('Original');

S1 = bowl_coords - focus_coords;
S2 = midpoint - (midpoint - [0, 0, 10]);
Theta = atan2(norm(cross(S1,S2)), dot(S1,S2));
t1_img_rot = imrotate3(t1_img, Theta, [-8, -2.5, -5]);
figure;
slice(double(t1_img_rot),size0(2)/2,size0(1)/2,size0(3)/2);
shading interp, colormap gray;
title('Rotated');