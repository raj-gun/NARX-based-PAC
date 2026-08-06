clear;clc;close all;

addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\NonSysID-i\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\Utils\');
%%
Fs = 1000; Ts = 1/Fs;
R=4;C=1;
%%
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);

%% =====================================================
%% PAC LF-sine HF-sine simple model
%% =====================================================

%%
tspan = 0:Ts:25;%(N*Ts-Ts);
N = length(tspan);
fftn = 4000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);

%% Nonsine PAC single LF coupling a single HF

am_lag = 16;
%----------
% s_LF = cos((2*pi*fL).*tspan + (f_phi(1)*pi/180));
%----------
LF_freq_1 = [6,7]; HF_freq_1 = [55,60];
n_smpls = 100; rng(200,"twister"); rng_seeds = randi([1,1e5],n_smpls,1); i = 78;
%----------
rng(rng_seeds(i), 'twister'); [s_LF, s_HF, ~, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_1);
rng(rng_seeds(i)+1000, 'twister'); [~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_1);
%=======================
% Non-sine PAC
m = 0.5; A = 200; C = 1*1e-6; 
[s_final, s_HF1_shft] = pac_general_1(s_LF, s_HF, A, C, m, am_lag);

% % Simple PAC
% m = 0.25; A = 50; 
% [s_final, s_HF1_shft] = pac_simple(s_LF, s_HF, A, m, am_lag);

% s_final = s_LF + s_HF1_shft;
%=======================

% pink_aug = (pink_aug_1 + pink_aug_2)./2;
disp(['LF = ', num2str(LF_freq_1), ', HF = ', num2str(HF_freq_1)]);


%%
% tm_frq_plt(s_final, Fs, fftn);
%% Visualise PAC signal
if am_lag~=0
    s_final = s_final(am_lag:end);
    tspan = tspan(am_lag:end);
    N = length(tspan);
    fftn = 1000;%Fs/N;
    w = 0:Fs/fftn:Fs-(Fs/fftn);
end
s_final_org = s_final;
s_final_org_fft = Ts.*fft(s_final_org, fftn);
% figure;subplot(2,1,1);plot(w, abs(s_final_org_fft));subplot(2,1,2);plot(w, angle(s_final_org_fft).*(180/pi));
% figure;plot(w, abs(s_final_org_fft));
%% Down sample

dwn_smpl_F = 250;
%s_final = fft_bndpss_flt( s_final, Fs, 0, dwn_smpl_F/2 );
s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
% s_final = lowpass_fir(s_final, (dwn_smpl_F/2)-2, Fs);
dwn_smpl = Fs/dwn_smpl_F;
s_final = s_final(1:dwn_smpl:N);
tspan = tspan(1:dwn_smpl:N);
Fs = dwn_smpl_F; Ts = 1/Fs;

pink = pink(1:dwn_smpl:N);

% N = length(tspan);
% fftn = 1000;%Fs/N;
% w = 0:Fs/fftn:Fs-(Fs/fftn);

%% Test single  sample of noise

%Trim the PAC signal to get a small segment
tm_windw = 10; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
s_final_trim = s_final(trim_ind)';%(600:1100)';

N = length(s_final_trim);
% figure;plot(tspan(trim_ind),s_final_trim);
%----------------------
fftn = length(s_final_trim); 
w = 0:Fs/fftn:Fs-(Fs/fftn);
tm_frq_plt(s_final_trim, Fs, fftn);
%% Add pink noise
pink = pink(trim_ind)';
s_final_trim = ( 3*( std(pink)/std(s_final_trim) ) ) .* s_final_trim;
sn_ratio = snr(s_final_trim,pink); disp(['SNR = ', num2str(sn_ratio), 'dB or ',  num2str(db2mag(sn_ratio))]);
%{1
s_final_trim = pink + s_final_trim;
tm_frq_plt(s_final_trim, Fs, fftn);
%}
%%

fL_vals = [2:1:15]; fH_vals = [20:1:100];

% fL_vals = [6:1:13]; fH_vals = [20:1:100];

% fL_vals = [3:1:10]; fH_vals = [HF_freq_1(1)-10:1:HF_freq_1(2)+10];


%% ------------------------ NARX based MISO PAC Comodulogram ------------------------------
filt_typ = {'sbp','sbp'}; % bw , sbp , guss
frq_bndw_LF = 1; 
frq_bndw_HF = 0.5;
%frq_bndw_HF = 0.1 & dwn_smpl_freq = 250; f<45 
%frq_bndw_HF = 0.2 & dwn_smpl_freq = 500; 45< f <65 
%frq_bndw_HF = 0.3 & dwn_smpl_freq = 500; xx< f <yy

disp(['frq_bndw_LF = ', num2str(frq_bndw_LF), ', frq_bndw_HF = ', num2str(frq_bndw_HF)]);

tic
[Comods , diff_comod, phs_data_mat, fL_grd, fH_grd, All_freq_comb_1, All_freq_comb, All_freq_comb_ARX_1, All_freq_comb_ARX_2, narx_pac_modls_1 , narx_pac_modls_2]...
    = pac_miso_Cmdg_mod_21(s_final_trim, fL_vals, fH_vals, Fs, 3, filt_typ, frq_bndw_LF, frq_bndw_HF);
toc

%%
fL_diff = mean(abs(diff(fL_vals))); fH_diff = mean(abs(diff(fH_vals))); 
rect_pos_box = @(LF_freq, HF_freq, fL_diff, fH_diff) [LF_freq(1)-fL_diff*0.5, HF_freq(1)-fH_diff*0.5, (abs(diff(LF_freq))*fL_diff)+1, (abs(diff(HF_freq))*1)+1]; 
rect_pos_1 = rect_pos_box(LF_freq_1, HF_freq_1, fL_diff, fH_diff);
% rect_pos_2 = rect_pos_box(LF_freq_1, HF_freq_2, fL_diff, fH_diff);
% rect_pos_3 = rect_pos_box(LF_freq_1, HF_freq_3, fL_diff, fH_diff);

figure; imagesc(fL_vals, fH_vals, Comods{1}); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
% rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
% rectangle('Position',rect_pos_3, 'EdgeColor','r', 'LineWidth', 1);

comod_D = diff_comod;
comod_D( comod_D > 0 ) = 1; comod_D( comod_D < 0 ) = -1; 
figure; imagesc(fL_vals, fH_vals, comod_D); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
% rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
% rectangle('Position',rect_pos_3, 'EdgeColor','r', 'LineWidth', 1);
%%  Post processing of results

%[HF_MI_score_1, HF_MI_score_final, Comod_intrmd] = SpuCup_intrmd(fL_vals, fH_vals, Comod, diff_comod);

% [Comod_harmonic_rmv, IF_harmonic_test_dat] = IF_harmonic_test (All_freq_comb_1, phs_data_mat, fL_vals, fH_vals, Comod, Ts);
% figure; imagesc(fL_vals, fH_vals, Comod_harmonic_rmv); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
% rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
% %rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
% % rectangle('Position',rect_pos_3, 'EdgeColor','r', 'LineWidth', 1);
% sgtitle('Commod after removing harmonics');

[HF_MI_score_1_2, HF_MI_score_final_2, Comod_intrmd] = SpuCup_intrmd_2(fL_vals, fH_vals, Comods{1}, diff_comod, All_freq_comb_1, phs_data_mat, Fs, 0);
figure; imagesc(fL_vals, fH_vals, Comod_intrmd); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
%rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
% rectangle('Position',rect_pos_3, 'EdgeColor','r', 'LineWidth', 1);
sgtitle('Commod after removing SC-i');

% [Comod_harmonic] = AM_var_commod(All_freq_comb_1, phs_data_mat, fL_vals, fH_vals);
% figure; imagesc(fL_vals, fH_vals, Comod_harmonic); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
% rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
% %rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);
% % rectangle('Position',rect_pos_3, 'EdgeColor','r', 'LineWidth', 1);
%% Local functions - PAC general

% Equation adapted from Jiang et al., (2015) NeuroImage
function [s_final, s_HF1_shft] = pac_general_1(s_LF, s_HF, a, c, m, delay_ind)
s_HF1 = m .* ( 1 - ( 1./(1 + exp(-a.*(s_LF-c))) ) ) .* s_HF;
if delay_ind ~= 0; s_HF1_shft = zeros(size(s_HF1)); s_HF1_shft(delay_ind:end) = s_HF1(1:end-delay_ind+1);
else; s_HF1_shft = s_HF1; end
s_final = s_LF + s_HF1_shft;
% s_final = s_final(delay_ind:3000+delay_ind-1); s_HF1_shft = s_HF1_shft(delay_ind:3000+delay_ind-1);
% figure;
% ax1=subplot(3,1,1);plot(s_HF1);hold on; plot(s_LF);
% ax2=subplot(3,1,2);plot(s_HF1_shft);hold on; plot(s_LF);
% ax3=subplot(3,1,3);plot(s_final); linkaxes([ax1,ax2,ax3],'x');
end
%-----------------------------------------------------

% Simplest form of PAC generaltion in electronics
%(J. Smith, Mathematics of the discrete Fourier transform (DFT). [North Charleston]: BookSurge, 2010.)
function [s_final, s_HF1_shft] = pac_simple(s_LF, s_HF, a, m, delay_ind)
s_HF1 = m .* (1 + a.*s_LF) .* s_HF;
if delay_ind~=0; s_HF1_shft = zeros(size(s_HF1)); s_HF1_shft(delay_ind:end) = s_HF1(1:end-delay_ind+1); else; s_HF1_shft = s_HF1; end
s_final = s_LF + s_HF1_shft;
% s_final = s_final(delay_ind:3000+delay_ind-1); s_HF1_shft = s_HF1_shft(delay_ind:3000+delay_ind-1);
% figure;
% ax1=subplot(3,1,1);plot(s_HF1);hold on; plot(s_LF);
% ax2=subplot(3,1,2);plot(s_HF1_shft);hold on; plot(s_LF);
% ax3=subplot(3,1,3);plot(s_final); linkaxes([ax1,ax2,ax3],'x');
%%
end



