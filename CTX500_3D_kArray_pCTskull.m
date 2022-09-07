
% =========================================================================
% CTX500_3D_kArray_pCTskull.m
%
% Run 3D k-wave acoustic simulations for the CTX-500 4-element transducer
% transducer modelled with kArray
% skull estimated from pCT generated from T1w-MRI
% example target: PCC
%
% simulation grid = 256x256x256
% grid spacing = 1x1x1 mm3
%
% script dependencies:
%   1) no1_prep.m
%   2) no2_runsim.m
%   3) fwhm2.m 
%   4) 
% Siti N. Yaakub | 25-Aug-2022 18:30:20
% =========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Change variables here to match your study requirements %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

% paths
base_dir = 'your_base_directory'; % your base directory
sim_output_dir = fullfile(base_dir, 'Sims'); % directory to save simulations to
% simulation pictures will be saved in the sim_output_dir

% subject-specific parameters
subj_id         = 'your_subject';       % subject ID
focus_coords    = [x1, y1, z1] + 1;     % point on the beam axis of the transducer [grid points]
bowl_coords     = [x2, y2, z2] + 1;     % centre of rear surface of transducer [grid points]
focus_depth     = 56;                   % to nearest mm

% ct and t1 filepaths (assumes these are saved in base_dir/Sims/t1 & base_dir/Sims/pct
ct_filename = 'path_to_ct/ct_filename.nii';
t1_filename = 'path_to_ct/t1_filename.nii';

% load correct version of NeuroFUS driving parameters
load driving_params_Isppa20.mat

% load parameters to use in your simulation 
% run savesimparams.m if you want to change any parameters
load simulation_params.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
base_dir = '/Users/sitiyaakub/Documents/Analysis/TUSMRS/';
output_dir = fullfile(base_dir, 'Sims');
load driving_params_Isppa20.mat

% subject-specific parameters
subj_id         = 'TS26_ERCA0870';      % subject ID
focus_coords    = [84, 91, 159] + 1;    % point on the beam axis of the transducer [grid points]
bowl_coords     = [66, 58, 227] + 1;     % centre of rear surface of transducer [grid points]
focus_depth     = 70;                   % to nearest mm

ct_filename = fullfile(base_dir, 'Sims', 'pct', ['sub-' subj_id '_T1w_n4_masked_pct.nii']);
t1_filename = fullfile(base_dir, 'Sims', 't1', ['sub-' subj_id '_T1w_n4_masked.nii']);


%% Prepare skull & check model
no1_prep.m 
% source, kgrid, model, t1_img, ppp, medium, rho_min,
% c_min, pulse_length, pulse_rep_freq, 
%% If 3D viewer of model looks OK, run simulation
no2_runsim.m 
% 

%% Summary
% --------------------
% SUMMARY
% --------------------
clear source;
disp(subj_id)
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
disp(['-6dB max width = ' num2str(max(field_focal_size_6dB)*1e3) ' mm'])
disp(' ')
disp('Copy to spreadsheet:')
disp([num2str(ppw) ';'  num2str(cfl) ';[' ...
    num2str(mx-shift_x) ',' num2str(my-shift_y) ',' num2str(mz-shift_z) '];' ...
    num2str(norm(bowl_coords-[mx,my,mz])*dx*1e3) ';' ...
    num2str(max_pressure * 1e-6) ';' num2str(MI) ';' num2str(Isppa) ';' ...
    num2str(p_focus * 1e-6) ';' num2str(isppa_focus) ';' ...
    num2str(focal_vol) ';' num2str(max(field_focal_size_6dB)*1e3)])

save(fullfile(sim_output_dir, [subj_id '_CTX500_3D_kArray_5Hz_rTUS_PCC.mat']));

% --------------------
% VISUALISATION
% --------------------

figure;
ax1 = axes;
imagesc(ax1, squeeze(t1_img(mx,:,:)), [50,500]);
% axis square; 
hold all;
ax2 = axes;
im2 = imagesc(ax2, squeeze(p(mx,:,:))*1e-6);
im2.AlphaData = 0.5;
% axis square;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; 
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude overlaid on T1')
saveas(gcf, fullfile(base_dir, 'Sims/Sim_pics', [subj_id '_5Hz_rTUS_PCC_sag.jpg']));

figure;
ax1 = axes;
imagesc(ax1, squeeze(t1_img(:,my,:)), [50,500]);
% axis square; 
hold all;
ax2 = axes;
im2 = imagesc(ax2, squeeze(p(:,my,:))*1e-6);
im2.AlphaData = 0.5;
% axis square;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; 
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude overlaid on T1')
saveas(gcf, fullfile(base_dir, 'Sims/Sim_pics', [subj_id '_5Hz_rTUS_PCC_cor.jpg']));

figure;
ax1 = axes;
imagesc(ax1, squeeze(t1_img(:,:,mz)), [50,500]);
% axis square; 
hold all;
ax2 = axes;
im2 = imagesc(ax2, squeeze(p(:,:,mz))*1e-6);
im2.AlphaData = 0.5;
% axis square;
linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; 
colormap(ax1,'gray')
colormap(ax2,'turbo')
set([ax1,ax2],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
xlabel(cb2, '[MPa]');
title(ax1,'Acoustic Pressure Amplitude overlaid on T1')
saveas(gcf, fullfile(base_dir, 'Sims/Sim_pics', [subj_id '_5Hz_rTUS_PCC_ax.jpg']));

%% for max p in brain (redo summary and plots if running this)
% find zz from plot of pressure in sagittal plane - x-axis coord below current max P
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

%% plots



%%
% % --------------------
% % THERMAL SIMULATION
% % --------------------
% 
% % convert the absorption coefficient to nepers/m
% alpha_np = db2neper(medium.alpha_coeff, medium.alpha_power) * ...
%     (2 * pi * freq).^medium.alpha_power;
% 
% % reshape the data, and calculate the volume rate of heat deposition
% Q = alpha_np .* p.^2 ./ (medium.density .* medium.sound_speed);
%
% run thermal simulation
% tic;
% % clear the input structures
% % clear medium source sensor;
% clear source;
% 
% % set the background temperature and heating term
% source2.Q = Q;
% source2.T0 = 37;
% 
% % define medium properties related to diffusion
% % ref: https://itis.swiss/virtual-population/tissue-properties/database
% % in skull
% medium2.density = rho_min + (rho_max - rho_min) * (model - 0) / (hu_max - 0);
% medium2.thermal_conductivity = zeros(size(model));
% medium2.thermal_conductivity(model > 0) = 0.32;
% medium2.specific_heat = 826 + (2524 - 826) * (1-((model - hu_min) / (hu_max - hu_min)));
% % in water
% medium2.density(model == 0)              = rho_min;     % [kg/m^3]
% medium2.thermal_conductivity(model == 0) = 0.6;      % [W/(m.K)]
% medium2.specific_heat(model == 0)        = 4178;     % [J/(kg.K)]
% 
% % create kWaveDiffusion object
% kdiff = kWaveDiffusion(kgrid, medium2, source2, [], 'DisplayUpdates', false, 'PlotSim', false);
% 
% % set source on time and off time
% on_time  = pulse_length;  % [s]
% off_time = 1/pulse_rep_freq - pulse_length;  % [s]
% num_bursts = stim_dur * pulse_rep_freq;
% % num_bursts = 2;
% % macaque protocol
% % pulse_length = 30e-3;
% % pulse_rep_freq = 10;
% % stim_dur = 40;
% % on_time  = pulse_length;  % [s]
% % off_time = 1/pulse_rep_freq - pulse_length;  % [s]
% % num_bursts = stim_dur * pulse_rep_freq;
% 
% % set time step size
% % dt = 0.1;
% dt = on_time/1;
% % dt = on_time/4;
% 
% maxT1 = zeros(size(model));
% for nb = 1:num_bursts
%     disp(['burst #' num2str(nb)])
%     % turn on heat source
%     kdiff.Q = Q;
%     % take time steps
%     kdiff.takeTimeStep(round(on_time / dt), dt);
%     
%     % store the current temperature field
%     T1 = kdiff.T;
%     
%     % turn off heat source and take time steps
%     kdiff.Q = 0;
%     kdiff.takeTimeStep(round(off_time / dt), dt);
%     
%     % store the current temperature field
%     T2 = kdiff.T;
% end
% toc;
% 
% disp(['Max temperature rise = ' num2str(max(T1(:))-source2.T0) ' degC'])
% save(fullfile(output_dir, [subj_id '_CTX500_3D_kArray_5Hz_rTUS_PCC_therm.mat']));
% 
% % --------------------
% % VISUALISATION
% % --------------------
% 
% % plot the thermal dose and lesion map (yz)
% figure;
% 
% % plot the acoustic pressure
% subplot(2, 2, 1);
% imagesc(kgrid.y_vec * 1e3, kgrid.z_vec * 1e3, squeeze(p(mx,:,:) * 1e-6).');
% h = colorbar;
% xlabel(h, '[MPa]');
% ylabel('z-position [mm]');
% xlabel('y-position [mm]');
% axis image;
% title('Acoustic Pressure Amplitude');
% 
% % plot the volume rate of heat deposition
% subplot(2, 2, 2);
% imagesc(kgrid.y_vec * 1e3, kgrid.z_vec * 1e3, squeeze(Q(mx,:,:) * 1e-7).');
% h = colorbar;
% xlabel(h, '[kW/cm^2]');
% ylabel('z-position [mm]');
% xlabel('y-position [mm]');
% axis image;
% title('Volume Rate Of Heat Deposition');
% 
% % plot the temperature after heating
% subplot(2, 2, 3);
% imagesc(kgrid.y_vec * 1e3, kgrid.z_vec * 1e3, squeeze(T1(mx,:,:)).');%, [37, 37.1]);
% h = colorbar;
% xlabel(h, '[degC]');
% ylabel('z-position [mm]');
% xlabel('y-position [mm]');
% axis image;
% title('Temperature After Heating');
% 
% % plot the temperature after cooling
% subplot(2, 2, 4);
% imagesc(kgrid.y_vec * 1e3, kgrid.z_vec * 1e3, squeeze(T2(mx,:,:)).');%, [37, 37.1]);
% h = colorbar;
% xlabel(h, '[degC]');
% ylabel('z-position [mm]');
% xlabel('y-position [mm]');
% axis image;
% title('Temperature After Cooling');
% 
% % set colormap and enlarge figure window
% % colormap(jet(256));
% colormap('turbo');
% scaleFig(2, 2);
% saveas(gcf, fullfile(base_dir, 'Sims/Sim_pics', [subj_id '_5Hz_rTUS_PCC_therm_yz.jpg']));
% 
% % plot the thermal dose and lesion map (xz)
% figure;
% 
% % plot the acoustic pressure
% subplot(2, 2, 1);
% imagesc(kgrid.x_vec * 1e3, kgrid.z_vec * 1e3, squeeze(p(:,my,:) * 1e-6).');
% h = colorbar;
% xlabel(h, '[MPa]');
% ylabel('z-position [mm]');
% xlabel('x-position [mm]');
% axis image;
% title('Acoustic Pressure Amplitude');
% 
% % plot the volume rate of heat deposition
% subplot(2, 2, 2);
% imagesc(kgrid.x_vec * 1e3, kgrid.z_vec * 1e3, squeeze(Q(:,my,:) * 1e-7).');
% h = colorbar;
% xlabel(h, '[kW/cm^2]');
% ylabel('z-position [mm]');
% xlabel('x-position [mm]');
% axis image;
% title('Volume Rate Of Heat Deposition');
% 
% % plot the temperature after heating
% subplot(2, 2, 3);
% imagesc(kgrid.x_vec * 1e3, kgrid.z_vec * 1e3, squeeze(T1(:,my,:)).');%, [37, 37.1]);
% h = colorbar;
% xlabel(h, '[degC]');
% ylabel('z-position [mm]');
% xlabel('x-position [mm]');
% axis image;
% title('Temperature After Heating');
% 
% % plot the temperature after cooling
% subplot(2, 2, 4);
% imagesc(kgrid.x_vec * 1e3, kgrid.z_vec * 1e3, squeeze(T2(:,my,:)).');%, [37, 37.1]);
% h = colorbar;
% xlabel(h, '[degC]');
% ylabel('z-position [mm]');
% xlabel('x-position [mm]');
% axis image;
% title('Temperature After Cooling');
% 
% % set colormap and enlarge figure window
% % colormap(jet(256));
% colormap('turbo');
% scaleFig(2, 2);
% saveas(gcf, fullfile(base_dir, 'Sims/Sim_pics', [subj_id '_5Hz_rTUS_PCC_therm_xz.jpg']));
% 
% % plot the thermal dose and lesion map (xy)
% figure;
% 
% % plot the acoustic pressure
% subplot(2, 2, 1);
% imagesc(kgrid.x_vec * 1e3, kgrid.y_vec * 1e3, (p(:,:,mz) * 1e-6).');
% h = colorbar;
% xlabel(h, '[MPa]');
% ylabel('y-position [mm]');
% xlabel('x-position [mm]');
% axis image;
% title('Acoustic Pressure Amplitude');
% 
% % plot the volume rate of heat deposition
% subplot(2, 2, 2);
% imagesc(kgrid.x_vec * 1e3, kgrid.y_vec * 1e3, (Q(:,:,mz) * 1e-7).');
% h = colorbar;
% xlabel(h, '[kW/cm^2]');
% ylabel('y-position [mm]');
% xlabel('x-position [mm]');
% axis image;
% title('Volume Rate Of Heat Deposition');
% 
% % plot the temperature after heating
% subplot(2, 2, 3);
% imagesc(kgrid.x_vec * 1e3, kgrid.y_vec * 1e3, T1(:,:,mz).');%, [37, 37.1]);
% h = colorbar;
% xlabel(h, '[degC]');
% ylabel('y-position [mm]');
% xlabel('x-position [mm]');
% axis image;
% title('Temperature After Heating');
% 
% % plot the temperature after cooling
% subplot(2, 2, 4);
% imagesc(kgrid.x_vec * 1e3, kgrid.y_vec * 1e3, T2(:,:,mz).');%, [37, 37.1]);
% h = colorbar;
% xlabel(h, '[degC]');
% ylabel('y-position [mm]');
% xlabel('x-position [mm]');
% axis image;
% title('Temperature After Cooling');
% 
% % set colormap and enlarge figure window
% % colormap(jet(256));
% colormap('turbo');
% scaleFig(2, 2);
% saveas(gcf, fullfile(base_dir, 'Sims/Sim_pics', [subj_id '_5Hz_rTUS_PCC_therm_xy.jpg']));
