clear all;clc;close all;

addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\Utils\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\Methods\Matlab_Code');
%%
Fs = 1000; Ts = 1/Fs;
R=4;C=1;
%%
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);
%% Load PAC data

load('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX-PAC paper\LFP\LFP data\LFP_HG_HFO.mat');
% s_final = lfpHG;
%or
s_final = lfpHFO;
N = length(s_final);
tspan = 0:Ts:(N*Ts-Ts);
fftn = 4000;%Fs/N;
tm_frq_plt(s_final, Fs, fftn); % Visualise PAC  in time-frequency plots
%% Down sample

% dwn_smpl_F = 500;
% %s_final = fft_bndpss_flt( s_final, Fs, 0, dwn_smpl_F/2 );
% s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
% % s_final = lowpass_fir(s_final, dwn_smpl_F/2, Fs);
% dwn_smpl = Fs/dwn_smpl_F;
% s_final = s_final(1:dwn_smpl:N);
% tspan = tspan(1:dwn_smpl:N);
% 
% Fs = dwn_smpl_F; Ts = 1/Fs;
% % N = length(tspan);
% % fftn = 1000;%Fs/N;
% % w = 0:Fs/fftn:Fs-(Fs/fftn);


%% Trim the PAC signal to get a small segment
tm_windw = 200; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
s_final_trim = s_final(trim_ind)';%(600:1100)';
N = length(s_final_trim);
figure;plot(tspan(trim_ind),s_final_trim);
%----------------------
fftn = length(s_final_trim); 
w = 0:Fs/fftn:Fs-(Fs/fftn);
% s_final_trim_fft = Ts.*fft(s_final_trim, fftn);figure;subplot(2,1,1);plot(w, abs(s_final_trim_fft));subplot(2,1,2);plot(w, angle(s_final_trim_fft).*(180/pi));
tm_frq_plt(s_final_trim, Fs, fftn);

fL_vals = [3:1:13]; fH_vals = [30:1:200];
% fL_vals = [3:0.5:20]; fH_vals = [30:1:200];
% fL_vals = [3:0.1:8]; fH_vals = [50:0.1:70];
% fL_vals = [3:0.1:3]; fH_vals = [67.4:0.1:67.4];


%% Evaluate PAC
[OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins] = modulationindex_directestimate_mod(s_final_trim', Fs,fL_vals,fH_vals,0.5,0.5,100);

phs_freq = reshape(Phs_Amp(flow_MI==8,:,:), size(Phs_Amp,2), size(Phs_Amp,3) );
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

% Save plot_data as HFO.mat