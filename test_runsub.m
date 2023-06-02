filepath = fullfile(pwd, '../');
t1_filename = fullfile(filepath, 'sub-test01_t1w.nii');
ct_filename = fullfile(filepath, 'sub-test01_pct.nii');
output_dir = filepath;
focus_coords_in = [99, 161, 202];  % point on the beam axis of the transducer [grid points]
bowl_coords_in = [90, 193, 262];   % centre of rear surface of transducer [grid points]
focus_depth = 60;               % to nearest mm

tussim_skull_3D(t1_filename, ct_filename, output_dir, ...
    focus_coords_in, bowl_coords_in, focus_depth)
tussim_skull_3D(t1_filename, ct_filename, output_dir, ...
    focus_coords_in, bowl_coords_in, focus_depth, ...
    'RunAcousticSim', true, 'RunThermalSim', true)

