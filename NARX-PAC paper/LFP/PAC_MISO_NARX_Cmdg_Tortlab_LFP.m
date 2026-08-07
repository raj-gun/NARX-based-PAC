%PAC_MISO_NARX_CMDG_TORTLAB_LFP Apply NARX-PAC to the Tortlab HFO LFP recording.
%   The script loads a 200 s HFO recording, evaluates coupling over 3-13 Hz
%   and 30-200 Hz, and plots the raw and post-processed comodulograms. Required
%   inputs are the LFP MAT-file, NonSysID-i, and the NARX-PAC utilities.
%
clear all;clc;close all;

addpath('\<path-to>\NonSysID-i\');
addpath('\<path-to>\NARX_PAC\');
addpath('\<path-to>\NARX_PAC\Utils\');
%% Set sampling and plotting parameters
Fs = 1000; Ts = 1/Fs;
R=4;C=1;
%% Define numerical helper functions
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);
%% Load the LFP recording

load('\<path-to>\NARX-PAC paper\LFP\LFP data\LFP_HG_HFO.mat');
%or
s_final = lfpHFO;
N = length(s_final);
tspan = 0:Ts:(N*Ts-Ts);
fftn = 4000;%Fs/N;
%% Downsample the LFP signal

dwn_smpl_F = 500;
s_final = lowpass_fir(s_final, (dwn_smpl_F/2)-2, Fs);
dwn_smpl = Fs/dwn_smpl_F;
s_final = s_final(1:dwn_smpl:N);
tspan = tspan(1:dwn_smpl:N);

Fs = dwn_smpl_F; Ts = 1/Fs;


%% Select the analysis segment
tm_windw = 200; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
s_final_trim = s_final(trim_ind)';%(600:1100)';
N = length(s_final_trim);
%----------------------
fftn = length(s_final_trim);
w = 0:Fs/fftn:Fs-(Fs/fftn);

fL_vals = [3:1:13]; fH_vals = [30:1:200];

%% Compute the NARX-based MISO PAC comodulogram
filt_typ = {'sbp','sbp'}; % bw , sbp , guss
frq_bndw_LF = 1; 
frq_bndw_HF = 0.5;

tic
[Comods , diff_comod, phs_data_mat, fL_grd, fH_grd, All_freq_comb_1, All_freq_comb, All_freq_comb_ARX_1, All_freq_comb_ARX_2, narx_pac_modls_1 , narx_pac_modls_2] = pac_miso_Cmdg_mod_21(s_final_trim, fL_vals, fH_vals, Fs, 3, filt_typ, frq_bndw_LF, frq_bndw_HF);
toc

% file_name = [num2str(Fs), 'Hz_', num2str(tm_windw), 's-', num2str(tm_itr), 'wndw_', num2str(frq_bndw_LF), '-', num2str(frq_bndw_HF), '-', filt_typ{1},'_', filt_typ{2}, '_LFcos'];
% file_dir = '/home/gunawardes/Documents/Matlab/CFC/Results/MISO_NARX/Tortlab/';
% save([file_dir , file_name, '.mat']);
%% Plot the raw and post-processed NARX-PAC maps
figure; imagesc(fL_vals, fH_vals, Comods{1}); colorbar; axis xy; set(gca, 'FontSize', 18);
figure; imagesc(fL_vals, fH_vals, Comods{2}); colorbar; axis xy; set(gca, 'FontSize', 18);
figure; imagesc(fL_vals, fH_vals, diff_comod); colorbar; axis xy; set(gca, 'FontSize', 18);

comod_D = diff_comod;
comod_D( comod_D > 0 ) = 1; comod_D( comod_D < 0 ) = -1; 
figure; imagesc(fL_vals, fH_vals, comod_D); colorbar; axis xy; set(gca, 'FontSize', 18);

[Comod_harmonic_rmv, IF_harmonic_test_dat] = IF_harmonic_test (All_freq_comb_1, phs_data_mat, fL_vals, fH_vals, Comods{1}, Ts);
figure; imagesc(fL_vals, fH_vals, Comod_harmonic_rmv); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
sgtitle('Commod after removing harmonics');

[HF_MI_score_1_2, HF_MI_score_final_2, Comod_intrmd] = SpuCup_intrmd_2(fL_vals, fH_vals, Comods{1}, diff_comod, All_freq_comb_1, phs_data_mat, Fs, 0);
figure; imagesc(fL_vals, fH_vals, Comod_intrmd); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
sgtitle('Commod after removing SC-i');


