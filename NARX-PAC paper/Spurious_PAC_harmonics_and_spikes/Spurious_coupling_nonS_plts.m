%SPURIOUS_COUPLING_NONS_PLTS Create the signal and method-comparison panels for Figure 20.
%   The script regenerates the noisy non-sinusoidal waveform, loads saved
%   NARX-PAC and conventional-method results, and formats the signal,
%   spectrum, and comodulograms for publication.
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
%% Configure the synthetic test signal
%% =====================================================

%% Configure the analysis interval
tspan = 0:Ts:25-Ts;%(N*Ts-Ts);
N = length(tspan);
fftn = 10000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);
%% Generate the non-sinusoidal waveform

LF = 7;
u = sin( (2*pi*LF).*tspan )./LF;
s_final = 10 ./ ( 1 + exp( -12.*( 0.5.*(1+5.*u) -0.7 ) )  );
am_lag = 0;


n_smpls = 100; rng(100,"twister"); rng_seeds = randi([1,1e5],n_smpls,1); 
rng( rng_seeds(60) ); [~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, [9.5,10.5], [55,65]);

rng(100,"twister");
SNR = 3;
s_final = ( SNR*( std(pink)/std(s_final) ) ) .* s_final;
sn_ratio = snr(s_final,pink); disp(['SNR = ', num2str(sn_ratio), 'dB or ',  num2str(db2mag(sn_ratio))]);
s_final = pink + s_final;


N = length(s_final);

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

OthrMthds_file_name = 'NonSine';
OthrMthds_file_dir = '\<path-to>\NARX-PAC paper\Spurious_PAC_harmonics_and_spikes\Data_Other_Methods\';
OthrMthds_plt_data = load([OthrMthds_file_dir, OthrMthds_file_name, '.mat']);
% plot_data = { {OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins} , ...
%               {GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM} };

MISO_NARX_file_dir = '\<path-to>\NARX-PAC paper\Spurious_PAC_harmonics_and_spikes\Data_NARX_PAC\';
MISO_NARX_file_name = 'NonSine';
%-----------------
MISO_NARX_plt_dat = load([MISO_NARX_file_dir,MISO_NARX_file_name,'.mat']);
MISO_NARX_plt_dat = { MISO_NARX_plt_dat.fL_vals, MISO_NARX_plt_dat.fH_vals, MISO_NARX_plt_dat.Comod_harmonic_rmv};

flow_MI = OthrMthds_plt_data.plot_data{1,1}{1,5};
fhigh_MI = OthrMthds_plt_data.plot_data{1,1}{1,6};

font_size = 24;
figure;
t3 = tiledlayout(1,3);
%---------
plots(w, s_final_fft, tspan, s_final, font_size, t3, 1, [1,1]);
%---------
Mthds_plt(OthrMthds_plt_data,MISO_NARX_plt_dat,font_size,t3, 2, [1,2], [1, 15, 15, 90], 84, -0.8);
%-------------------------


%% =====================================================
%% Local functions
%% =====================================================


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

function Mthds_plt(OthrMthds_plt_data,MISO_NARX_plt_dat,font_size,tile_tag, tile_no, tile_spn, axis_lim, tt_pos, ylab_pos)
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
penny_comod = OthrMthds_plt_data.plot_data{1,2}{1,2}; penny_comod = log(penny_comod);
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