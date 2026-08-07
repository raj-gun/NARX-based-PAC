clear;clc;close all;
addpath('\<path-to>\NARX_PAC\Utils');
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
%% Generate PAC signal

am_lag = 16;
%----------
% s_LF = cos((2*pi*fL).*tspan + (f_phi(1)*pi/180));
%----------
LF_freq = [6,7]; HF_freq = [55,60];
n_smpls = 100; rng(200,"twister"); rng_seeds = randi([1,1e5],n_smpls,1); i = 15;
%----------
rng(rng_seeds(i), 'twister'); [s_LF, s_HF, ~, ~] = pink_noise_LF_HF(N, Fs, LF_freq, HF_freq);
rng(rng_seeds(i)+1000, 'twister'); [~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, LF_freq, HF_freq);
%=======================
% Non-sine PAC
m = 0.5; A = 200; C = 1*1e-6; am_lag = 16;
[s_final, s_HF1_shft] = pac_general_1(s_LF, s_HF, A, C, m, am_lag);

%=======================

% pink_aug = (pink_aug_1 + pink_aug_2)./2;
disp(['LF = ', num2str(LF_freq), ', HF = ', num2str(HF_freq)]);


s_final = ( 3*( std(pink)/std(s_final) ) ) .* s_final;
sn_ratio = snr(s_final,pink); disp(['SNR = ', num2str(sn_ratio), 'dB or ',  num2str(db2mag(sn_ratio))]);
s_final = pink + s_final;


%%

tm_frq_plt(s_final, Fs, fftn);
%% Visualise PAC signal
if am_lag~=0
    s_final = s_final(am_lag:end);
    tspan = tspan(am_lag:end);
    N = length(tspan);
    %fftn = 1000;%Fs/N;
    w = 0:Fs/fftn:Fs-(Fs/fftn);
end
s_final_fft = Ts.*fft(s_final, fftn);
figure;subplot(2,1,1);plot(w, abs(s_final_fft));subplot(2,1,2);plot(w, angle(s_final_fft).*(180/pi));
% figure;plot(w, abs(s_final_org_fft));
%%

OthrMthds_file_name = 'PAC_OthrMthds_6-7_55-60_SNR_3';
OthrMthds_file_dir = '\<path-to>\NARX-PAC paper\6-7_PAC_55-65\';
OthrMthds_plt_data = load([OthrMthds_file_dir, OthrMthds_file_name, '.mat']);
% plot_data = { {OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins} , ...
%               {GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM} };

MISO_NARX_file_dir = '\<path-to>\NARX-PAC paper\6-7_PAC_55-65\';
MISO_NARX_file_name = '6-7_55-65_pinknoise_SNR_3';
%-----------------
MISO_NARX_plt_dat = load([MISO_NARX_file_dir,MISO_NARX_file_name,'.mat']);
MISO_NARX_plt_dat = { MISO_NARX_plt_dat.fL_vals, MISO_NARX_plt_dat.fH_vals, MISO_NARX_plt_dat.Comod, MISO_NARX_plt_dat.Comod_intrmd, MISO_NARX_plt_dat.comod_D};

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

%% PAC general

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

%% Spike train

function [spike_train] = spike_signal(mean_interval,N,jitter,amplitude,width_samples,Fs)
%----------------- Generate pink noise -----------------------
white = randn(1, N);
f = fft(white);
frequencies = [ 0:Fs/N:(Fs/2)-(Fs/N) , fliplr(0:Fs/N:(Fs/2)) ];
scaling = 1 ./ sqrt(frequencies); % 1/f amplitude decay
scaling(1)=0; scaling(end)=0;
f = f .* scaling;
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
imagesc( MISO_NARX_plt_dat{1} , MISO_NARX_plt_dat{2} , MISO_NARX_plt_dat{4} );
xlabel('Low Frequency (Hz)');
ttl_4 = title({'NARX-based PAC';' '}); ttl_4.Position(2)=ttl_3.Position(2);
axis([3, 15, 30, 95]);
colorbar; axis xy; set(gca, 'FontSize', font_size);
%----------
%----------
% nexttile(t2)
% imagesc( MISO_NARX_plt_dat{1} , MISO_NARX_plt_dat{2} , MISO_NARX_plt_dat{4} );
% ttl_5 = title({'Penny et. al. 2008';' '}); ttl_5.Position(2)=89;
% xlabel('Low Frequency (Hz)');
% %ylh_2 = ylabel('High Frequency (Hz)'); ylh_2.Position(1)=ylh_1.Position(1);
% axis([3, 15, 30, 95]);
% colorbar; axis xy; set(gca, 'FontSize', font_size);
% %----------
% nexttile(t2)
% imagesc( MISO_NARX_plt_dat{1} , MISO_NARX_plt_dat{2} , MISO_NARX_plt_dat{5} );
% xlabel('Low Frequency (Hz)');
% ttl_4 = title({'NARX-based PAC';' '}); ttl_4.Position(2)=ttl_3.Position(2);
% axis([3, 15, 30, 95]);
% colorbar; axis xy; set(gca, 'FontSize', font_size);
end