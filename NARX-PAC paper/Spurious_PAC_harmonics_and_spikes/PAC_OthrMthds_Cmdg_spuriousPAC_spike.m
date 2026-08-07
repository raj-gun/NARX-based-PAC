clear all;clc;close all;

addpath('/home/gunawardes/Documents/Matlab/CFC/NARX_PAC/Utils/');
addpath('/home/gunawardes/Documents/Matlab/CFC/Methods/Matlab_Code/');
%%
Fs = 1000; Ts = 1/Fs;
R=4;C=1;
%%
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);

%% =====================================================
%% Spurious PAC
%% =====================================================

%%
tspan = 0:Ts:25-Ts;%(N*Ts-Ts);
N = length(tspan);
fftn = 4000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);

%% Gaussian-shaped Spike train

LF = 5;
amplitude = 3;            % In units of std dev of background
width_ms = 10;            % FWHM in milliseconds
width_samples = round((width_ms / 1000) * Fs);
% Spike train parameters
mean_interval = 1000/LF;      % Mean interval in ms (e.g., 100 ms = 10 Hz)
jitter = 20;              % Jitter in ms
rng(600);
spike_train = spike_signal(mean_interval,N,jitter,amplitude,width_samples,Fs);
s_final = spike_train;
tm_frq_plt( s_final, Fs, N);

rng(570);
[~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, [9.5,10.5], [55,65]);
s_final = 1.19.*s_final;
sn_ratio = snr(s_final,pink); disp(['SNR = ', num2str(sn_ratio)]);
s_final = pink + s_final;

%% Down sample

% dwn_smpl_F = 250;
% %s_final = fft_bndpss_flt( s_final, Fs, 0, dwn_smpl_F/2 );
% s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
% % s_final = lowpass_fir(s_final, dwn_smpl_F/2, Fs);
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
tm_windw = 20; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
s_final_trim = s_final(trim_ind)';%(600:1100)';

N = length(s_final_trim);
figure;plot(tspan(trim_ind),s_final_trim);
%----------------------
fftn = length(s_final_trim); 
w = 0:Fs/fftn:Fs-(Fs/fftn);
% s_final_trim_fft = Ts.*fft(s_final_trim, fftn);figure;subplot(2,1,1);plot(w, abs(s_final_trim_fft));subplot(2,1,2);plot(w, angle(s_final_trim_fft).*(180/pi));
tm_frq_plt(s_final_trim, Fs, fftn);

fL_vals = [1:1:20]; fH_vals = [15:1:90];
% fL_vals = [1:0.5:15]; fH_vals = [18:1:50];
% fL_vals = [5:0.5:20]; fH_vals = [30:1:100];


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
% Use the following code snippet to save the workspace for plotting

% save_file_name = 'spike_train';
% save_file_dir = '\<path-to>\NARX-PAC paper\Spurious_PAC_harmonics_and_spikes\Data_Other_Methods\';
% save([save_file_dir, save_file_name], 'plot_data');
% ===============================
%%

%% =====================================================
%% Local functions
%% =====================================================S

%% Spike train
function [spike_train] = spike_signal(mean_interval,N,jitter,amplitude,width_samples,Fs)
%----------------- Generate pink noise -----------------------
white = randn(1, N);
f = fft(white);
frequencies = [ 0:Fs/N:(Fs/2)-(Fs/N) , fliplr(0:Fs/N:(Fs/2)) ];
scaling = 1 ./ sqrt(frequencies); % 1/f amplitude decay
scaling(1)=0; scaling(end)=0;
f = f .* scaling(1:end-1);
pink_noise = real(ifft(f));
%-------------------------------------------------------------
% pink_noise = dsp.ColoredNoise('Color','pink','SamplesPerFrame',N,'NumChannels',1);
% bg_signal = pink_noise();
bg_signal = pink_noise;
sigma = width_samples / (2*sqrt(2*log(2))); % Convert FWHM to standard deviation for Gaussian
intervals_ms = mean_interval + (rand(1, ceil(N/(Fs*mean_interval/1000))) - 0.5)*2*jitter;
intervals_samples = round(intervals_ms / 1000 * Fs);
spike_times = cumsum(intervals_samples);
spike_times(spike_times > N) = [];

% Generate spike train
spike_train = zeros(1, N);
for i = 1:length(spike_times)
    center = spike_times(i);
    % Create Gaussian window centered at this spike time
    window = -round(3*sigma):round(3*sigma);
    gauss = exp(-0.5 * (window / sigma).^2);
    idx = center + window;
    valid = idx > 0 & idx <= N;
    spike_train(idx(valid)) = spike_train(idx(valid)) + amplitude * std(bg_signal) * gauss(valid);
end

end