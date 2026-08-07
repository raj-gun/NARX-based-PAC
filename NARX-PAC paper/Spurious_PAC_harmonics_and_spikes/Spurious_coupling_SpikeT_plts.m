clear;clc;close all;
addpath('\<path-to>\NARX_PAC\Utils\');
%%
Fs = 1000; Ts = 1/Fs;
%%
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);

%% =====================================================
%% PAC LF-sine HF-sine simple model
%% =====================================================

%%
tspan = 0:Ts:25-Ts;%(N*Ts-Ts);
N = length(tspan);
fftn = 10000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);
%% Non-sinusoidal signal

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

%%
tm_frq_plt(s_final, Fs, fftn);

s_final_fft = Ts.*fft(s_final, fftn);
figure;subplot(2,1,1);plot(w, abs(s_final_fft));subplot(2,1,2);plot(w, angle(s_final_fft).*(180/pi));
% figure;plot(w, abs(s_final_org_fft));
%%

OthrMthds_file_name = 'spike_train';
OthrMthds_file_dir = '\<path-to>\NARX-PAC paper\Spurious_PAC_harmonics_and_spikes\Data_Other_Methods\';
OthrMthds_plt_data = load([OthrMthds_file_dir, OthrMthds_file_name, '.mat']);
% plot_data = { {OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins} , ...
%               {GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM} };

MISO_NARX_file_dir = '\<path-to>\NARX-PAC paper\Spurious_PAC_harmonics_and_spikes\Data_NARX_PAC\';
MISO_NARX_file_name = 'spike_train';
%-----------------
MISO_NARX_plt_dat = load([MISO_NARX_file_dir,MISO_NARX_file_name,'.mat']);
MISO_NARX_plt_dat = { MISO_NARX_plt_dat.fL_vals, MISO_NARX_plt_dat.fH_vals, MISO_NARX_plt_dat.Comod};

flow_MI = OthrMthds_plt_data.plot_data{1,1}{1,5};
fhigh_MI = OthrMthds_plt_data.plot_data{1,1}{1,6};

font_size = 24;
figure;
t3 = tiledlayout(1,3);
%---------
plots(w, s_final_fft, tspan, s_final, font_size, t3, 1, [1,1]);
%---------
Mthds_plt(OthrMthds_plt_data,MISO_NARX_plt_dat,font_size,t3, 2, [1,2], [1, 20, 15, 90], 84, -1.2);
%-------------------------


%% =====================================================
%% Local functions
%% =====================================================

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

%% Plotting
function plots(w, s_final_fft, tspan, s_final, font_size, tile_tag, tile_no, tile_spn)
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
axis([4.5 4.5+5 min(s_final) max(s_final)]); xticks([4.5, 4.5+5]); xticklabels({'0','5'});
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

function Mthds_plt(OthrMthds_plt_data,MISO_NARX_plt_dat,font_size,tile_tag, tile_no, tile_spn, axis_lim, tt_pos, ylab_pos)
flow_MI = OthrMthds_plt_data.plot_data{1,1}{1,5};
fhigh_MI = OthrMthds_plt_data.plot_data{1,1}{1,6};
t2 = tiledlayout(tile_tag,2,2, 'TileSpacing','loose', 'Padding', 'loose');
t2.Layout.Tile = tile_no;
t2.Layout.TileSpan = tile_spn;
%----------
nexttile(t2)
imagesc(flow_MI, fhigh_MI, OthrMthds_plt_data.plot_data{1,1}{1,1}');
ttl_1 = title({'Ozkurt et. al. 2011';' '}); ttl_1.Position(2)=tt_pos;
ylh_1 = ylabel('High Frequency (Hz)'); ylh_1.Position(1)=ylab_pos;
axis(axis_lim);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
nexttile(t2)
canolty_MI_comod = OthrMthds_plt_data.plot_data{1,1}{1,2}; canolty_MI_comod( canolty_MI_comod > 1e-1 ) = 0;
imagesc(flow_MI, fhigh_MI, canolty_MI_comod');
ttl_2 = title({'Canolty et. al. 2010';' '}); ttl_2.Position(2)=ttl_1.Position(2);
axis(axis_lim);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
nexttile(t2)
penny_comod = OthrMthds_plt_data.plot_data{1,2}{1,2}; penny_comod = log10(penny_comod);
imagesc(OthrMthds_plt_data.plot_data{1,2}{1,4}, OthrMthds_plt_data.plot_data{1,2}{1,5}, penny_comod');
ttl_3 = title({'Penny et. al. 2008';' '}); ttl_3.Position(2)=tt_pos;
xlabel('Low Frequency (Hz)');
ylh_2 = ylabel('High Frequency (Hz)'); ylh_2.Position(1)=ylh_1.Position(1);
axis(axis_lim);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
nexttile(t2)
imagesc( MISO_NARX_plt_dat{1} , MISO_NARX_plt_dat{2} , MISO_NARX_plt_dat{3} );
xlabel('Low Frequency (Hz)');
ttl_4 = title({'NARX-based PAC';' '}); ttl_4.Position(2)=ttl_3.Position(2);
axis(axis_lim);
colorbar; axis xy; set(gca, 'FontSize', font_size);
end