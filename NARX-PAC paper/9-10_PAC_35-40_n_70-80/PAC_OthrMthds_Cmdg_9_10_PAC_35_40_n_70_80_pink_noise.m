clear all;clc;close all;

addpath('\<path-to>\NARX_PAC\Utils\');
addpath('\<path-to>\Methods\Matlab_Code');
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
%% Nonsine PAC single LF coupling with two distinct HFs

am_lag_1 = 16;
am_lag_2 = 90;
% am_lag_3 = 20;
%----------
% s_LF = cos((2*pi*fL).*tspan + (f_phi(1)*pi/180));
%----------
LF_freq_1 = [9,10]; 
% HF_freq_1 = [60,75]; HF_freq_2 = [40,50]; HF_freq_3 = [20,30]+10; i=2;
% HF_freq_1 = [65,75]+5; HF_freq_2 = [40,45]-5; HF_freq_3 = [20,30]+10; i = 4;%

HF_freq_1 = [65,75]+5; HF_freq_2 = [40,45]-5; i = 4; % Working
% HF_freq_1 = [65,75]; HF_freq_2 = [40,45]; HF_freq_3 = [20,30]+10; i=2; %Not working

% HF_freq_1 = [60,75]+10; HF_freq_2 = [40,50]-10; HF_freq_3 = [20,30]+10; i = 4;
% HF_freq_1 = [65,75]; HF_freq_2 = [40,45]; HF_freq_3 = [20,30]+10; i = 5;

n_smpls = 100; rng(100,"twister"); rng_seeds = randi([1,1e5],n_smpls,1); 
%----------
rng(rng_seeds(i), 'twister'); [s_LF_1, s_HF_1, pink, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_1);
[~, s_HF_2, ~, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_2);
% [~, s_HF_3, ~, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_3);
rng(rng_seeds(i)+1000, 'twister'); [~, ~, pink1, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_2);
%=======================
% Non-sine PAC
% m = 0.75; A = 100; C = 1*1e-6;
% [~, s_HF1_shft_1] = pac_general_1(s_LF_1, s_HF_1, A, C, m, am_lag_1);
% [~, s_HF1_shft_2] = pac_general_1(s_LF_1, s_HF_2, A, C, m, am_lag_2);
% s_final = s_LF_1 +  s_HF1_shft_1 + s_HF1_shft_2;

% m1 = 0.5; A1 = 50;
% m2 = 0.5; A2 = 50;
% [~, s_HF1_shft_1] = pac_simple(s_LF_1, s_HF_1, A1, m1, am_lag_1);
% [~, s_HF1_shft_2] = pac_simple(s_LF_1, s_HF_2, A2, m2, am_lag_2);
% s_final = s_LF_1 +  s_HF1_shft_1 + s_HF1_shft_2;

m1 = 0.5; A1 = 50;
m2 = 0.25; A2 = 150;
m3 = 0.5; A3 = 50;
[~, s_HF1_shft_1] = pac_simple(s_LF_1, s_HF_1, A1, m1, am_lag_1);
[~, s_HF1_shft_2] = pac_simple(s_LF_1, s_HF_2, A2, m2, am_lag_2);
% [~, s_HF1_shft_3] = pac_simple(s_LF_1, s_HF_3, A3, m3, am_lag_3);
s_final = s_LF_1 +  s_HF1_shft_1 + s_HF1_shft_2;% + s_HF1_shft_3;
%=======================
am_lag = max([am_lag_1,am_lag_2]);

% pink_aug = (pink_aug_1 + pink_aug_2)./2;
disp(['LF = ', num2str(LF_freq_1), ', HF1 = ', num2str(HF_freq_1), ', HF2 = ', num2str(HF_freq_2)]);

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

% dwn_smpl_F = 250;
% %s_final = fft_bndpss_flt( s_final, Fs, 0, dwn_smpl_F/2 );
% s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
% % s_final = lowpass_fir(s_final, (dwn_smpl_F/2)-2, Fs);
% dwn_smpl = Fs/dwn_smpl_F;
% s_final = s_final(1:dwn_smpl:N);
% tspan = tspan(1:dwn_smpl:N);
% Fs = dwn_smpl_F; Ts = 1/Fs;
% 
% pink = pink(1:dwn_smpl:N);

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


%% Evaluate PAC
[OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins] = modulationindex_directestimate_mod(s_final_trim', Fs,fL_vals,fH_vals,0.5,0.5,100);

phs_freq = reshape(Phs_Amp(flow_MI==6,:,:), size(Phs_Amp,2), size(Phs_Amp,3) );
% phs_freq = phs_freq./max(phs_freq')';
figure;imagesc(phs_bins, fhigh_MI, phs_freq ); axis xy;

figure
imagesc(flow_MI, fhigh_MI, OzktMI');  colorbar;
axis xy
set(gca, 'FontSize', 18);
% xlabel('Phase Frequency (Hz)');  ylabel('Amplitude Frequency (Hz)');
% title('Canolty Modulation index')

figure
imagesc(flow_MI, fhigh_MI, CanltyMI');  colorbar;
axis xy
set(gca, 'FontSize', 18);
% xlabel('Phase Frequency (Hz)');  ylabel('Amplitude Frequency (Hz)');
% title('Direct estimator')

figure
imagesc(flow_MI, fhigh_MI, TortMI');  colorbar;
axis xy
set(gca, 'FontSize', 18);
% xlabel('Phase Frequency (Hz)');  ylabel('Amplitude Frequency (Hz)');
% title('Direct estimator')

[GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM] = general_linear_index_mod(s_final_trim', Fs,fL_vals,fH_vals,0.5,0.5);

figure
imagesc(flow_GLM, fhigh_GLM, GLM_org');  colorbar;
axis xy
set(gca, 'FontSize', 18);
% xlabel('Phase Frequency (Hz)');  ylabel('Envelope Frequency (Hz)');


figure;
imagesc(flow_GLM, fhigh_GLM, GLM_2');  colorbar;
axis xy
set(gca, 'FontSize', 18);
% xlabel('Phase Frequency (Hz)');  ylabel('Envelope Frequency (Hz)');


figure
imagesc(flow_GLM, fhigh_GLM, GLM_robust');  colorbar;
axis xy
set(gca, 'FontSize', 18);
% xlabel('Phase Frequency (Hz)');  ylabel('Envelope Frequency (Hz)');

plot_data = { {OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins} , ...
              {GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM} };

% ===============================
% Use the following code snippet to save the workspace for plotting using
% Synthetic_exmpl_9_10__35_40_n_70_80_SNR_3_plts.m

% save_file_name = 'PAC_OthrMthds_9_10_PAC_35_40_n_70_80';
% save_file_dir = '\<path-to>\NARX-PAC paper\9-10_PAC_35-40_n_70-80\';
% save([save_file_dir, save_file_name], 'plot_data');
% ===============================
%%
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



