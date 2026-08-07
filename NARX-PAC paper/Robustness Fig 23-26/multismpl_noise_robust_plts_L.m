%MULTISMPL_NOISE_ROBUST_PLTS_L Create the noise-robustness comparison in Figure 25.
%   The script loads saved median results for the two-band experiment at SNR 2
%   and SNR 1, compares NARX-PAC with the conventional methods, and marks the
%   known coupled frequency bands.
%
clear all;clc;close all;

%% Define coupling-region annotations
rect_pos_box = @(LF_freq, HF_freq, fL_diff, fH_diff) [LF_freq(1)-fL_diff*0.5, HF_freq(1)-fH_diff*0.5, (abs(diff(LF_freq))*1)+1, (abs(diff(HF_freq))*1)+1]; 

%% Load the SNR 2 results
MISO_NARX_dir = '\<path-to>\NARX-PAC paper\9-10_PAC_35-40_n_70-80\Median_100_samples\Plot data\';
load([MISO_NARX_dir,'F1_2_10s_wrk_PP.mat']);
Comod_NARX_PAC = Comods_mat;
Comod_D_NARX_PAC = Comods_D_mat;
Comods_PP_intrmd = Comods_intrmd_mat; 

Othr_Mthds_dir = '\<path-to>\NARX-PAC paper\9-10_PAC_35-40_n_70-80\Median_100_samples\Plot data\';
load([Othr_Mthds_dir,'PAC_OthrMthds_9-10__35-40__70-80_SNR_2.mat']);

%% Plot the SNR 2 comparison
fL_diff = mean(abs(diff(fL_vals))); fH_diff = mean(abs(diff(fH_vals))); 
rect_pos_1 = rect_pos_box(LF_freq_1, HF_freq_1, fL_diff, fH_diff);
rect_pos_2 = rect_pos_box(LF_freq_1, HF_freq_2, fL_diff, fH_diff);

figure; tiledlayout(2,4, 'TileSpacing','loose', 'Padding', 'loose');
font_size = 19;
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comods_mat_OzktMI,3)'); caxis([0,3e-3]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
ylabel('High Frequency (Hz)');
title({'Ozkurt et. al. 2011';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comods_mat_CanltyMVL,3)'); caxis([0,0.04]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
title({'Canolty et. al. 2010';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comods_mat_TortMI,3)'); caxis([0,0.025]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
title({'Tort et. al. 2008';' '});
%
nexttile;
imagesc(flow_GLM, fhigh_GLM, median(Comods_mat_GLM_robust,3)'); caxis([0,3e-3]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
xlabel('Low Frequency (Hz)');
title({'Penny et. al. 2008';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comod_NARX_PAC,3)); caxis([0,1]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
ylabel('High Frequency (Hz)');
xlabel('Low Frequency (Hz)');
title({'NARX-based PAC';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comods_PP_intrmd,3)); caxis([0,1]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
xlabel('Low Frequency (Hz)');
title({'NARX-based PAC PP';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comod_D_NARX_PAC,3)); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
xlabel('Low Frequency (Hz)');
title({'NARX-based PAC D'; ' '});

%% Load the SNR 1 results
MISO_NARX_dir = '\<path-to>\NARX-PAC paper\9-10_PAC_35-40_n_70-80\Median_100_samples\Plot data\';
load([MISO_NARX_dir,'F1_1_10s_wrk_PP.mat']);
Comod_NARX_PAC = Comods_mat;
Comod_D_NARX_PAC = Comods_D_mat;
Comods_PP_intrmd = Comods_intrmd_mat; 

Othr_Mthds_dir = '\<path-to>\NARX-PAC paper\9-10_PAC_35-40_n_70-80\Median_100_samples\Plot data\';
load([Othr_Mthds_dir,'PAC_OthrMthds_9-10__35-40__70-80_SNR_1.mat']);

%% Plot the SNR 1 comparison
fL_diff = mean(abs(diff(fL_vals))); fH_diff = mean(abs(diff(fH_vals))); 
rect_pos_1 = rect_pos_box(LF_freq_1, HF_freq_1, fL_diff, fH_diff);
rect_pos_2 = rect_pos_box(LF_freq_1, HF_freq_2, fL_diff, fH_diff);

figure; tiledlayout(2,4, 'TileSpacing','loose', 'Padding', 'loose');
font_size = 19;
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comods_mat_OzktMI,3)'); caxis([0,3e-3]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
ylabel('High Frequency (Hz)');
title({'Ozkurt et. al. 2011';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comods_mat_CanltyMVL,3)'); caxis([0,0.04]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
title({'Canolty et. al. 2010';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comods_mat_TortMI,3)'); caxis([0,0.025]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
title({'Tort et. al. 2008';' '});
%
nexttile;
imagesc(flow_GLM, fhigh_GLM, median(Comods_mat_GLM_robust,3)'); caxis([0,3e-3]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
xlabel('Low Frequency (Hz)');
title({'Penny et. al. 2008';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comod_NARX_PAC,3)); caxis([0,1]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
ylabel('High Frequency (Hz)');
xlabel('Low Frequency (Hz)');
title({'NARX-based PAC';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comods_PP_intrmd,3)); caxis([0,1]); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
xlabel('Low Frequency (Hz)');
title({'NARX-based PAC PP';' '});
%
nexttile;
imagesc(fL_vals, fH_vals, median(Comod_D_NARX_PAC,3)); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
xlabel('Low Frequency (Hz)');
title({'NARX-based PAC D'; ' '});
