clear all;clc;close all;

%% Load saved data
MISO_NARX_dir = 'C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX-PAC paper\6-7_PAC_55-65\Median_100_samples\Plot data\';
dat_1_2 = load([MISO_NARX_dir,'F2_3_2s_wrk_PP.mat']); 
Comod_NARX_PAC_1_2 = {dat_1_2.Comods_intrmd_mat, dat_1_2.fH_vals, dat_1_2.fL_vals};

dat_1_3 = load([MISO_NARX_dir,'F2_3_3s_wrk_PP.mat']);
Comod_NARX_PAC_1_3 = {dat_1_3.Comods_intrmd_mat, dat_1_3.fH_vals, dat_1_3.fL_vals};

dat_1_5 = load([MISO_NARX_dir,'F2_3_5s_wrk_PP.mat']);
Comod_NARX_PAC_1_5 = {dat_1_5.Comods_intrmd_mat, dat_1_5.fH_vals, dat_1_5.fL_vals};

LF_freq = dat_1_2.LF_freq_1; HF_freq = dat_1_2.HF_freq_1;    

fL_diff = mean(abs(diff(dat_1_2.fL_vals))); fH_diff = mean(abs(diff(dat_1_2.fH_vals)));

clear dat_1_2 dat_1_3 dat_1_5;

MISO_NARX_dir = 'C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX-PAC paper\9-10_PAC_35-40_n_70-80\Median_100_samples\Plot data\';
dat_2_2 = load([MISO_NARX_dir,'F1_3_2s_wrk_PP.mat']);
Comod_NARX_PAC_2_2 = {dat_2_2.Comods_intrmd_mat, dat_2_2.fH_vals, dat_2_2.fL_vals};

dat_2_3 = load([MISO_NARX_dir,'F1_3_3s_wrk_PP.mat']);
Comod_NARX_PAC_2_3 = {dat_2_3.Comods_intrmd_mat, dat_2_3.fH_vals, dat_2_3.fL_vals};

dat_2_5 = load([MISO_NARX_dir,'F1_3_5s_wrk_PP.mat']);
Comod_NARX_PAC_2_5 = {dat_2_5.Comods_intrmd_mat, dat_2_5.fH_vals, dat_2_5.fL_vals};

LF_freq_1 = dat_2_2.LF_freq_1; HF_freq_1 = dat_2_2.HF_freq_1; HF_freq_2 = dat_2_2.HF_freq_2;    

fL_diff_1 = mean(abs(diff(dat_2_2.fL_vals))); fH_diff_1 = mean(abs(diff(dat_2_2.fH_vals)));


clear dat_2_2 dat_2_3 dat_2_5;

%% Plots
 
rect_pos_box = @(LF_freq, HF_freq, fL_diff, fH_diff) [LF_freq(1)-fL_diff*0.5, HF_freq(1)-fH_diff*0.5, (abs(diff(LF_freq))*1)+1, (abs(diff(HF_freq))*1)+1]; 

rect_pos = rect_pos_box(LF_freq, HF_freq, fL_diff, fH_diff);

rect_pos_1 = rect_pos_box(LF_freq_1, HF_freq_1, fL_diff_1, fH_diff_1);
rect_pos_2 = rect_pos_box(LF_freq_1, HF_freq_2, fL_diff_1, fH_diff_1);

%===========================================

figure; tiledlayout(2,3, 'TileSpacing', 'loose', 'Padding', 'loose');
font_size = 19;

% Compute the plotted matrices first
Z_1_2 = median(Comod_NARX_PAC_1_2{1,1}, 3);
Z_1_3 = median(Comod_NARX_PAC_1_3{1,1}, 3);
Z_1_5 = median(Comod_NARX_PAC_1_5{1,1}, 3);

Z_2_2 = median(Comod_NARX_PAC_2_2{1,1}, 3);
Z_2_3 = median(Comod_NARX_PAC_2_3{1,1}, 3);
Z_2_5 = median(Comod_NARX_PAC_2_5{1,1}, 3);

% Common colour range for all plots
allZ = [Z_1_2(:); Z_1_3(:); Z_1_5(:); Z_2_2(:); Z_2_3(:); Z_2_5(:)];
allZ = allZ(~isnan(allZ));   % ignore NaNs if present

commonCLim = [0, 1];

%
nexttile;
imagesc(Comod_NARX_PAC_1_2{1,3}, Comod_NARX_PAC_1_2{1,2}, median(Comod_NARX_PAC_1_2{1,1},3) ); caxis(commonCLim); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos, 'EdgeColor','r', 'LineWidth', 1);
% ylabel({'\bf{One FO}' ; ' ' ; '\rm{High Frequency (Hz)}'});
ylabel('High Frequency (Hz)');
%xlabel('Low Frequency (Hz)');
title({'Signal length of 2 seconds'});
% set(gca,'Units','normalized'); titleHandle = get( gca ,'Title' ); pos  = get( titleHandle , 'position' ); pos1 = [pos(1) pos(2)-7 pos(3)]; set( titleHandle , 'position' , pos1 );
%
nexttile;
imagesc(Comod_NARX_PAC_1_3{1,3}, Comod_NARX_PAC_1_3{1,2}, median(Comod_NARX_PAC_1_3{1,1},3) ); caxis(commonCLim); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos, 'EdgeColor','r', 'LineWidth', 1);
% ylabel('High Frequency (Hz)');
%xlabel('Low Frequency (Hz)');
title({'Signal length of 3 seconds'});
% set(gca,'Units','normalized'); titleHandle = get( gca ,'Title' ); pos  = get( titleHandle , 'position' ); pos1 = [pos(1) pos(2)-7 pos(3)]; set( titleHandle , 'position' , pos1 );
%
nexttile;
imagesc(Comod_NARX_PAC_1_5{1,3}, Comod_NARX_PAC_1_5{1,2}, median(Comod_NARX_PAC_1_5{1,1},3) ); caxis(commonCLim); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos, 'EdgeColor','r', 'LineWidth', 1);
% ylabel('High Frequency (Hz)');
%xlabel('Low Frequency (Hz)');
title({'Signal length of 5 seconds'});
% set(gca,'Units','normalized'); titleHandle = get( gca ,'Title' ); pos  = get( titleHandle , 'position' ); pos1 = [pos(1) pos(2)-7 pos(3)]; set( titleHandle , 'position' , pos1 );
%
nexttile;
imagesc(Comod_NARX_PAC_2_2{1,3}, Comod_NARX_PAC_2_2{1,2}, median(Comod_NARX_PAC_2_2{1,1},3) ); caxis(commonCLim); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
% ylabel({'\bf{Two FOs}' ; ' ' ; '\rm{High Frequency (Hz)}'});
ylabel('High Frequency (Hz)');
xlabel('Low Frequency (Hz)');
%title({'5 sec.';' '});
% set(gca,'Units','normalized'); titleHandle = get( gca ,'Title' ); pos  = get( titleHandle , 'position' ); pos1 = [pos(1) pos(2)-7 pos(3)]; set( titleHandle , 'position' , pos1 );
%
nexttile;
imagesc(Comod_NARX_PAC_2_3{1,3}, Comod_NARX_PAC_2_3{1,2}, median(Comod_NARX_PAC_2_3{1,1},3) ); caxis(commonCLim); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
% ylabel('High Frequency (Hz)');
xlabel('Low Frequency (Hz)');
%title({'5 sec.';' '});
% set(gca,'Units','normalized'); titleHandle = get( gca ,'Title' ); pos  = get( titleHandle , 'position' ); pos1 = [pos(1) pos(2)-7 pos(3)]; set( titleHandle , 'position' , pos1 );
%
nexttile;
imagesc(Comod_NARX_PAC_2_5{1,3}, Comod_NARX_PAC_2_5{1,2}, median(Comod_NARX_PAC_2_5{1,1},3) ); caxis(commonCLim); colorbar; axis xy; set(gca, 'FontSize', font_size); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
% ylabel('High Frequency (Hz)');
xlabel('Low Frequency (Hz)');
%title({'5 sec.';' '});
% set(gca,'Units','normalized'); titleHandle = get( gca ,'Title' ); pos  = get( titleHandle , 'position' ); pos1 = [pos(1) pos(2)-7 pos(3)]; set( titleHandle , 'position' , pos1 );

