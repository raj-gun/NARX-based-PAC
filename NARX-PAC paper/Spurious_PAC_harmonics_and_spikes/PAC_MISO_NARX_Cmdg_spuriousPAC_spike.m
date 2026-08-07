clear all;clc;close all;

addpath('\<path-to>\NonSysID-i\');
addpath('\<path-to>\NARX_PAC\');
addpath('\<path-to>\NARX_PAC\Utils\');
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

dwn_smpl_F = 250;
s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
dwn_smpl = Fs/dwn_smpl_F;
s_final = s_final(1:dwn_smpl:N);
tspan = tspan(1:dwn_smpl:N);
Fs = dwn_smpl_F; Ts = 1/Fs;

pink = pink(1:dwn_smpl:N);


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
tm_frq_plt(s_final_trim, Fs, fftn);

fL_vals = [1:1:20]; fH_vals = [15:1:90];


%% ------------------------ NARX based MISO PAC Comodulogram ------------------------------
filt_typ = {'sbp','sbp'}; % bw , sbp , guss
frq_bndw_LF = 0.5;
frq_bndw_HF = 0.5;
disp(['frq_bndw_LF = ', num2str(frq_bndw_LF), ', frq_bndw_HF = ', num2str(frq_bndw_HF)]);

tic
[Comods , diff_comod, phs_data_mat, fL_grd, fH_grd, All_freq_comb_1, All_freq_comb, All_freq_comb_ARX_1, All_freq_comb_ARX_2, narx_pac_modls_1 , narx_pac_modls_2]...
    = pac_miso_Cmdg_mod_21(s_final_trim, fL_vals, fH_vals, Fs, 3, filt_typ, frq_bndw_LF, frq_bndw_HF);
toc
%% Post processing of results


[Comod_harmonic_rmv, IF_harmonic_test_dat] = IF_harmonic_test(All_freq_comb_1, phs_data_mat, fL_vals, fH_vals, Comods{1}, Ts);
figure; imagesc(fL_vals, fH_vals, Comod_harmonic_rmv); colorbar; axis xy; set(gca, 'FontSize', 18);
sgtitle('Commod after removing harmonics');


%%
figure; imagesc(fL_vals, fH_vals, Comods{1}); colorbar; axis xy; set(gca, 'FontSize', 18);



%%
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