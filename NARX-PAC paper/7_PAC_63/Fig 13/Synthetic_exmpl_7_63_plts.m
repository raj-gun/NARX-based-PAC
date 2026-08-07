%SYNTHETIC_EXMPL_7_63_PLTS Assemble the Figure 13 synthetic-signal and method-comparison panels.
%   The script regenerates the noisy 7 Hz by 63 Hz PAC signal, loads saved
%   outputs from NARX-PAC and the conventional methods, and formats the
%   signal, spectrum, comodulograms, and phase-amplitude summaries.
%
clear;clc;close all;
addpath('\<path-to>\NARX_PAC\Utils');
%% Set sampling and plotting parameters
Fs = 1000; Ts = 1/Fs;
%% Define numerical helper functions
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);

%% =====================================================
%% PAC LF-sine HF-sine simple model
%% =====================================================

%% Configure the analysis interval
tspan = 0:Ts:25-Ts;%(N*Ts-Ts);
N = length(tspan);
fftn = 10000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);
%% Generate PAC signal

fL=7; fH=63;
f_phi = [39.4212  290.5146];
am_lag = 35;
%=======================
s_LF = cos((2*pi*fL).*tspan + (f_phi(1)*pi/180));
s_HF = cos((2*pi*fH).*tspan + (f_phi(2)*pi/180));

LF_freq = [6.95,7.05]; HF_freq = [62.95,63.05];
rng(130);%rng(130);
[~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, LF_freq, HF_freq);
%=======================
% Non-sine PAC
m = 1;
[s_final, s_HF1_shft] = pac_general_1(s_LF, s_HF, 6, 1*1e-6, m, am_lag);
%=======================
s_final = ( 3*( std(pink)/std(s_final) ) ) .* s_final;
s_final = pink + s_final;


disp([fL,fH]); disp([f_phi,am_lag]);
%% Visualise the current results

tm_frq_plt(s_final, Fs, fftn);
%% Visualise PAC signal
if am_lag~=0
    s_final = s_final(am_lag:end);
    tspan = tspan(am_lag:end);
    N = length(tspan);
    w = 0:Fs/fftn:Fs-(Fs/fftn);
end
s_final_fft = Ts.*fft(s_final, fftn);
figure;subplot(2,1,1);plot(w, abs(s_final_fft));subplot(2,1,2);plot(w, angle(s_final_fft).*(180/pi));

%% Load saved method-comparison results

OthrMthds_file_name = 'PAC_OthrMthds_7_63';
OthrMthds_file_dir = '\<path-to>\NARX-PAC paper\7_PAC_63\Fig 13\';
OthrMthds_plt_data = load([OthrMthds_file_dir, OthrMthds_file_name, '.mat']);
% plot_data = { {OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins} , ...
%               {GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM} };

MISO_NARX_file_dir = '\<path-to>\NARX-PAC paper\7_PAC_63\Fig 13\';
MISO_NARX_file_name = '7-63_pinknoise';
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

%% Phase angle
%Trim the PAC signal to get a small segment
tm_windw = 10; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
s_final_trim = s_final(trim_ind)';%(600:1100)';
std_cos = sqrt(0.5);
s_final_trim = ( (s_final_trim-mean(s_final_trim)).* ( std_cos/std(s_final_trim) ) );

if am_lag~=0
    s_LF = s_LF(am_lag:end);
    s_HF1_shft = s_HF1_shft(am_lag:end);
end
%----------------------------------------
LF_phs_mod = angle(hilbert(s_LF));
HF_env_mod = abs(hilbert(s_HF1_shft));
nbins = 100;
[bin_centers, amp_means] = pac_hist(nbins , LF_phs_mod , HF_env_mod , true);

%----------------------------------------

str = '#ff0000b2'; color_lf = [sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 , 0.7]; % low frequency signal color
str = '#0000ffb2'; color_hf = [sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 , 0.7]; % high freqeuncy singal color
figure;ax=subplot(1,1,1);plot(s_LF(1:length(trim_ind)), 'Color', color_lf,'LineWidth',2);axis([0,inf,-inf,inf]);hold on;
plot(s_HF1_shft(1:length(trim_ind)), 'Color', color_hf,'LineWidth',2);axis([0,inf,-inf,inf]);
box(ax,'off');set(get(ax, 'YAxis'), 'Visible', 'off');set(get(ax, 'XAxis'), 'Visible', 'off'); %axis([0,400,-inf,inf]);

figure;plot(trim_ind,s_final_trim);

%% =====================================================
%% Local functions
%% =====================================================

function [bin_centers, amp_means] = pac_hist(nbins,LF_phs,HF_env, plt)
%PAC_HIST Average the fast amplitude within slow-phase bins.
%   NBINS sets the number of phase bins, LF_PHS contains slow phase in radians,
%   HF_ENV is the fast-amplitude envelope, and PLT enables an optional bar
%   plot. Outputs are BIN_CENTERS and the corresponding AMP_MEANS.
edges = linspace(-pi, pi, nbins+1); % Phase bins
amp_means = zeros(1, nbins);        % To store mean amplitudes
% Bin the phase data and compute mean amplitude in each bin
for k = 1:nbins
    indices = LF_phs >= edges(k) & LF_phs < edges(k+1);
    amp_means(k) = mean(HF_env(indices));
end
% For visualization: bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
if plt
    figure;
    bar(rad2deg(bin_centers), amp_means, 'FaceColor', [0.2 0.6 0.8], 'EdgeAlpha', 0);
    set(gca,"TickLabelInterpreter",'latex'); set(gca, 'XTick', [-180 -90 0 90 180]); set(gca, 'XTickLabel', {'$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'});
    set(gca,'YTick',[]);
    set(gca, 'FontSize', 20);
end
end

%% PAC general

% Equation adapted from Jiang et al., (2015) NeuroImage
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

% Simplest form of PAC generaltion in electronics
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
%   at half maximum, and FS is the sampling rate. SPIKE_TRAIN is the resulting
%   noisy periodic-transient signal.
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

%% Plotting
function plots(w, s_final_fft, tspan, s_final, font_size, tile_tag, tile_no, tile_spn)
%PLOTS Draw the publication time-trace and magnitude-spectrum panel.
%   W and S_FINAL_FFT provide the frequency axis and spectrum; TSPAN and
%   S_FINAL provide the time-domain data. FONT_SIZE and the TILE_* arguments
%   control placement in the parent tiled layout. This helper returns no data.
%---------------
box_top = max(abs(s_final_fft))*1.0;
%---------------
str = '#006801e3'; color_raw = [sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 , 1]; % raw singal color
str = '#D95319'; color_fft = [sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 , 1]; % Magnitude spectrum line color
%---------------
t1=tiledlayout(tile_tag,4,1);
t1.Layout.Tile = tile_no;
t1.Layout.TileSpan = tile_spn;
%---------------
ax1 = nexttile(t1);
plot(tspan, s_final, 'Color', color_raw,'LineWidth',2); % Raw signal plot
axis([4.5 5.75 min(s_final) max(s_final)]); xticks([4.5, 5, 5.5]); xticklabels({'0','0.5','1'});
xlabel('Time (sec)');
ylabel('$s(t)$', 'Interpreter','latex', 'Rotation', 0);
%---
ax4 = nexttile(t1,2,[3,1]);  % Magntitude spectrum
plot(w, abs(s_final_fft), 'Color', color_fft, 'LineWidth',2);
axis([0 125 0 box_top]);
xlabel('Frequency (Hz)');
ylabel('$\left| S(\omega) \right|$', 'Interpreter','latex',  'Rotation', 0);
%---------------
set(ax1,'YTick',[]); set(ax1, 'FontSize', font_size); box(ax1,'off');
set(ax4,'YTick',[]); set(ax4, 'FontSize', font_size); box(ax4,'off');
end

function Mthds_plt(OthrMthds_plt_data,MISO_NARX_plt_dat,font_size,tile_tag, tile_no, tile_spn)
%MTHDS_PLT Draw conventional and NARX-PAC comodulogram panels.
%   OTHRMTHDS_PLT_DATA and MISO_NARX_PLT_DAT contain saved method outputs;
%   FONT_SIZE and TILE_* control the tiled layout. Any remaining arguments
%   set scenario-specific axis limits and label positions. This helper
%   produces plots and returns no data.
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
axis([4, 10, 30, 95]);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
nexttile(t2)
imagesc(flow_MI, fhigh_MI, OthrMthds_plt_data.plot_data{1,1}{1,2}');
ttl_2 = title({'Canolty et. al. 2010';' '}); ttl_2.Position(2)=ttl_1.Position(2);
axis([4, 10, 30, 95]);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
nexttile(t2)
imagesc(OthrMthds_plt_data.plot_data{1,2}{1,4}, OthrMthds_plt_data.plot_data{1,2}{1,5}, OthrMthds_plt_data.plot_data{1,2}{1,2}');
ttl_3 = title({'Penny et. al. 2008';' '}); ttl_3.Position(2)=89;
xlabel('Low Frequency (Hz)');
ylh_2 = ylabel('High Frequency (Hz)'); ylh_2.Position(1)=ylh_1.Position(1);
axis([4, 10, 30, 95]);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
nexttile(t2)
imagesc( MISO_NARX_plt_dat{1} , MISO_NARX_plt_dat{2} , MISO_NARX_plt_dat{3} );
xlabel('Low Frequency (Hz)');
ttl_4 = title({'NARX-based PAC';' '}); ttl_4.Position(2)=ttl_3.Position(2);
axis([4, 10, 30, 95]);
colorbar; axis xy; set(gca, 'FontSize', font_size);
end