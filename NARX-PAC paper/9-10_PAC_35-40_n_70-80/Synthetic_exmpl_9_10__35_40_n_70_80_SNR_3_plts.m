%SYNTHETIC_EXMPL_9_10__35_40_N_70_80_SNR_3_PLTS Create the signal and method-comparison panels for Figure 16.
%   The script regenerates the synthetic PAC signal containing two coupled
%   fast-frequency bands, loads saved NARX-PAC and conventional-method results,
%   and formats the time trace, spectrum, and comodulograms for publication.
%
clear;clc;close all;
addpath('\<path-to>\NARX_PAC\Utils\');
%% Set sampling and plotting parameters
Fs = 1000; Ts = 1/Fs;
%% Define numerical helper functions
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);

%% =====================================================
%% Configure the synthetic PAC experiment
%% =====================================================

%% Configure the analysis interval
tspan = 0:Ts:25-Ts;%(N*Ts-Ts);
N = length(tspan);
fftn = 10000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);
%% Generate PAC signal

am_lag_1 = 16;
am_lag_2 = 90;
%----------
%----------
LF_freq_1 = [6,7]; HF_freq_1 = [65,75]+5; HF_freq_2 = [40,45]-5;

n_smpls = 100; rng(100,"twister"); rng_seeds = randi([1,1e5],n_smpls,1); i = 1;

rng(rng_seeds(i), 'twister'); [s_LF_1, s_HF_1, pink, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_1);
[~, s_HF_2, ~, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_2);
rng(rng_seeds(i)+1000, 'twister'); [~, ~, pink1, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_2);
%=======================
% Non-sine PAC
m1 = 0.5; A1 = 50;
m2 = 0.25; A2 = 150;
m3 = 0.5; A3 = 50;
[~, s_HF1_shft_1] = pac_simple(s_LF_1, s_HF_1, A1, m1, am_lag_1);
[~, s_HF1_shft_2] = pac_simple(s_LF_1, s_HF_2, A2, m2, am_lag_2);
s_final = s_LF_1 +  s_HF1_shft_1 + s_HF1_shft_2;
%=======================
am_lag = max(am_lag_1,am_lag_2);

disp(['LF = ', num2str(LF_freq_1), ', HF1 = ', num2str(HF_freq_1), ', HF2 = ', num2str(HF_freq_2) ]);


s_final = ( 3*( std(pink)/std(s_final) ) ) .* s_final;
sn_ratio = snr(s_final,pink); disp(['SNR = ', num2str(sn_ratio), 'dB or ',  num2str(db2mag(sn_ratio))]);
s_final = pink + s_final;

%% Inspect the signal in the time and frequency domains

tm_frq_plt(s_final, Fs, fftn);
%% Align the signal and compute its spectrum
if am_lag~=0
    s_final = s_final(am_lag:end);
    tspan = tspan(am_lag:end);
    N = length(tspan);
    w = 0:Fs/fftn:Fs-(Fs/fftn);
end
s_final_fft = Ts.*fft(s_final, fftn);
figure;subplot(2,1,1);plot(w, abs(s_final_fft));subplot(2,1,2);plot(w, angle(s_final_fft).*(180/pi));
%% Load saved method-comparison results

OthrMthds_file_name = 'PAC_OthrMthds_9-10__35-40__70-80_SNR_3';
OthrMthds_file_dir = '\<path-to>\NARX-PAC paper\9-10_PAC_35-40_n_70-80\';
OthrMthds_plt_data = load([OthrMthds_file_dir, OthrMthds_file_name, '.mat']);
% plot_data = { {OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins} , ...
%               {GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM} };

MISO_NARX_file_dir = '\<path-to>\NARX-PAC paper\9-10_PAC_35-40_n_70-80\';
MISO_NARX_file_name = '9-10_35-40_70-80_pinknoise_SNR_3';
%-----------------
MISO_NARX_plt_dat = load([MISO_NARX_file_dir,MISO_NARX_file_name,'.mat']);
MISO_NARX_plt_dat = { MISO_NARX_plt_dat.fL_vals, MISO_NARX_plt_dat.fH_vals, MISO_NARX_plt_dat.Comod_intrmd};

flow_MI = OthrMthds_plt_data.plot_data{1,1}{1,5};
fhigh_MI = OthrMthds_plt_data.plot_data{1,1}{1,6};

font_size = 24;
figure;
t3 = tiledlayout(1,3);
%---------
plots(w, s_final_fft, tspan, s_final, font_size, t3, 1, [1,1]);
%---------
Mthds_plt(OthrMthds_plt_data,MISO_NARX_plt_dat,font_size,t3, 2, [1,2]);
%-------------------------


%% =====================================================
%% Local functions
%% =====================================================

%% PAC signal-generation functions

% Equation adapted from Jiang et al. (2015), NeuroImage.
function [s_final, s_HF1_shft] = pac_general_1(s_LF, s_HF, a, c, m, delay_ind)
%PAC_GENERAL_1 Generate PAC using non-sinusoidal amplitude modulation.
%   S_LF and S_HF are the slow and fast components; A and C control the
%   logistic modulation shape, M scales the fast component, and DELAY_IND
%   delays the modulated fast component in samples. Outputs are the composite
%   signal S_FINAL and the shifted amplitude-modulated component S_HF1_SHFT.
s_HF1 = m .* ( 1 - ( 1./(1 + exp(-a.*(s_LF-c))) ) ) .* s_HF;
if delay_ind ~= 0; s_HF1_shft = zeros(size(s_HF1)); s_HF1_shft(delay_ind:end) = s_HF1(1:end-delay_ind+1);
else; s_HF1_shft = s_HF1; end
s_final = s_LF + s_HF1_shft;
end
%-----------------------------------------------------

% Basic linear amplitude-modulation model of PAC.
%(J. Smith, Mathematics of the discrete Fourier transform (DFT). [North Charleston]: BookSurge, 2010.)
function [s_final, s_HF1_shft] = pac_simple(s_LF, s_HF, a, m, delay_ind)
%PAC_SIMPLE Generate PAC using linear amplitude modulation.
%   S_LF and S_HF are the slow and fast components; A is the modulation
%   depth, M scales the fast component, and DELAY_IND applies a sample delay.
%   Outputs are the composite signal S_FINAL and the shifted modulated fast
%   component S_HF1_SHFT.
s_HF1 = m .* (1 + a.*s_LF) .* s_HF;
if delay_ind~=0; s_HF1_shft = zeros(size(s_HF1)); s_HF1_shft(delay_ind:end) = s_HF1(1:end-delay_ind+1); else; s_HF1_shft = s_HF1; end
s_final = s_LF + s_HF1_shft;
end

%% Spike train

function [spike_train] = spike_signal(mean_interval,N,jitter,amplitude,width_samples,Fs)
%SPIKE_SIGNAL Generate a jittered Gaussian spike train on pink noise.
%   MEAN_INTERVAL and JITTER are in milliseconds, N is the sample count,
%   AMPLITUDE sets the spike height, WIDTH_SAMPLES is the Gaussian full width
%   at half maximum, and FS is the sampling rate. SPIKE_TRAIN contains the
%   pink-noise background and Gaussian spikes.
%----------------- Generate pink noise -----------------------
white = randn(1, N);
f = fft(white);
frequencies = [ 0:Fs/N:(Fs/2)-(Fs/N) , fliplr(0:Fs/N:(Fs/2)) ];
scaling = 1 ./ sqrt(frequencies); % 1/f amplitude decay
scaling(1)=0; scaling(end)=0;
f = f .* scaling;
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

%% Plotting functions
function plots(w, s_final_fft, tspan, s_final, font_size, tile_tag, tile_no, tile_spn)
%PLOTS Plot the time-domain signal and its magnitude spectrum.
%   W and S_FINAL_FFT provide the frequency axis and spectrum; TSPAN and
%   S_FINAL provide the time-domain data. FONT_SIZE sets the text size, while
%   TILE_TAG, TILE_NO, and TILE_SPN define the position in the parent tiled
%   layout. The function has no output arguments.
%---------------
box_top = max(abs(s_final_fft))*1.0;
%---------------
str = '#006801e3'; color_raw = [sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 , 1]; % Raw-signal colour
str = '#D95319'; color_fft = [sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 , 1]; % Magnitude-spectrum line colour
%---------------
t1=tiledlayout(tile_tag,4,1);
t1.Layout.Tile = tile_no;
t1.Layout.TileSpan = tile_spn;
%---------------
ax1 = nexttile(t1);
plot(tspan, s_final, 'Color', color_raw,'LineWidth',2); % Raw signal plot
axis([4.5 5.75 min(s_final) max(s_final)]); xticks([4.5, 5, 5.5]); xticklabels({'0','0.5','1'});
xlabel('Time (sec)');
ylabel('$z(t)$', 'Interpreter','latex', 'Rotation', 0);
%---
ax4 = nexttile(t1,2,[3,1]);  % Magntitude spectrum
plot(w, abs(s_final_fft), 'Color', color_fft, 'LineWidth',2);
axis([0 125 0 box_top]);
xlabel('Frequency (Hz)');
ylabel('$\left| Z(\omega) \right|$', 'Interpreter','latex',  'Rotation', 0);
%---------------
set(ax1,'YTick',[]); set(ax1, 'FontSize', font_size); box(ax1,'off');
set(ax4,'YTick',[]); set(ax4, 'FontSize', font_size); box(ax4,'off');
end

function Mthds_plt(OthrMthds_plt_data,MISO_NARX_plt_dat,font_size,tile_tag, tile_no, tile_spn)
%MTHDS_PLT Plot conventional and NARX-PAC comodulograms.
%   OTHRMTHDS_PLT_DATA and MISO_NARX_PLT_DAT contain the saved method results.
%   FONT_SIZE sets the text size; TILE_TAG, TILE_NO, and TILE_SPN define the
%   tiled-layout position. AXIS_LIM sets the displayed frequency ranges, while
%   TT_POS and YLAB_POS adjust the title and y-axis label positions. The
%   function has no output arguments.
flow_MI = OthrMthds_plt_data.plot_data{1,1}{1,5};
fhigh_MI = OthrMthds_plt_data.plot_data{1,1}{1,6};
t2 = tiledlayout(tile_tag,2,2, 'TileSpacing','loose', 'Padding', 'loose');
t2.Layout.Tile = tile_no;
t2.Layout.TileSpan = tile_spn;
%----------
nexttile(t2)
imagesc(flow_MI, fhigh_MI, OthrMthds_plt_data.plot_data{1,1}{1,1}');
ttl_1 = title({'Ozkurt et. al. 2011';' '}); ttl_1.Position(2)=89;
ylh_1 = ylabel('High Frequency (Hz)'); ylh_1.Position(1)=1.4;
axis([3, 15, 30, 95]);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
nexttile(t2)
imagesc(flow_MI, fhigh_MI, OthrMthds_plt_data.plot_data{1,1}{1,2}');
ttl_2 = title({'Canolty et. al. 2010';' '}); ttl_2.Position(2)=ttl_1.Position(2);
axis([3, 15, 30, 95]);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
nexttile(t2)
imagesc(OthrMthds_plt_data.plot_data{1,2}{1,4}, OthrMthds_plt_data.plot_data{1,2}{1,5}, OthrMthds_plt_data.plot_data{1,2}{1,2}');
ttl_3 = title({'Penny et. al. 2008';' '}); ttl_3.Position(2)=89;
xlabel('Low Frequency (Hz)');
ylh_2 = ylabel('High Frequency (Hz)'); ylh_2.Position(1)=ylh_1.Position(1);
axis([3, 15, 30, 95]);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
nexttile(t2)
imagesc( MISO_NARX_plt_dat{1} , MISO_NARX_plt_dat{2} , MISO_NARX_plt_dat{3} );
xlabel('Low Frequency (Hz)');
ttl_4 = title({'NARX-based PAC';' '}); ttl_4.Position(2)=ttl_3.Position(2);
axis([3, 15, 30, 95]);
colorbar; axis xy; set(gca, 'FontSize', font_size);
end