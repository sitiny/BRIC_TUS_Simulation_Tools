subj_id         = 'SSTO0540';      % subject ID
focus_coords    = [79, 139, 150]*2 + 1;    % point on the beam axis of the transducer [grid points]
bowl_coords     = [14, 170, 194]*2 + 1;    % centre of rear surface of transducer [grid points]
focus_depth     = 75;                   % to nearest mm
transducer = 'CTX500';

base_dir = ['~/sub-' subj_id];
output_dir = fullfile(base_dir, 'plan-nacc_ppw6');
ct_filename = fullfile(base_dir, ['sub-' subj_id '_pCT_0.5mm.nii.gz']);
t1_filename = fullfile(base_dir, ['sub-' subj_id '_T1w_0.5mm.nii.gz']);

pressure = 97360;
phases = [0	275.3	190.6	105.8];

%tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth, transducer);

tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth, transducer, 'RunAcousticSim', true, 'SourcePressure', pressure, 'SourcePhase', phases);
