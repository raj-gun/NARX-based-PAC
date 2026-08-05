
%SIMULATE_MISO_NARX_PAC Inspect one identified NARX-PAC model in detail.
%   This post-identification analysis script loads a saved synthetic-data
%   result, selects one detected frequency pair, simulates the identified
%   model, decomposes it into the canonical low-frequency and amplitude-
%   modulated components, evaluates preferred phase, and inspects the
%   sideband-to-carrier ratios used by the modulation index (equation (15)).
%   Required loaded variables include the frequency grids and comodulogram,
%   All_freq_comb, narx_pac_modls_1, Fs, Ts, N, frq_bndw_LF,
%   frq_bndw_HF, filt_typ, s_final_trim, and the synthetic reference
%   components used in the comparison plots.
%   This file is a script, not a function. It does not perform the complete
%   PAC grid identification. Update the absolute paths and selected pair index
%   before use. NonSysID-i is required for model_simulation_i and related
%   model utilities.

%% Reset the workspace and configure external paths
clear all;clc;close all;
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\NonSysID-i\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\Utils');
%% Load a saved NARX-PAC experiment and inspect its comodulogram
file_dir = 'C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\Results\MISO_NARX\SyntheticData\';
file_name = '7-63_pinknoise';
load([file_dir,file_name]);
figure; imagesc(fL_vals, fH_vals, Comod); colorbar; axis xy; set(gca, 'FontSize', 18); %axis([3,20, 35,95,-inf,inf]);
figure; surf(fL_grd, fH_grd, Comod); shading interp; %axis([2,20, 25.5,90,-inf,inf]);
%% Define numerical helpers used by the analysis
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);
freq_pos = @(freq,Fs,fftn) floor(freq*fftn/Fs) + 1;
%% Select a detected frequency pair and display its identified terms
j=2;
pos_freq_comp = All_freq_comb(j,[1,2])
model = narx_pac_modls_1{j,1};
Model_terms = model{1,16}.Row;
Model_coefficients = model{1,16}.theta;
Model_tbl = table(Model_terms(1:end-1,1), Model_coefficients(1:end-1))
phi = [0,0];
fftn_dsrd = (Fs/0.1);
tspan_nofrf = 0:Ts:(fftn_dsrd+1000)*Ts;%-6+Ts:Ts:6-Ts;%tspan_sysid;%
u_fft = [ cos( 2.*pi.*pos_freq_comp(1).*tspan_nofrf + phi(1) )' , cos( 2.*pi.*pos_freq_comp(2).*tspan_nofrf + phi(2) )' ];

fftn = 25000;
w = 0:Fs/fftn:Fs-(Fs/fftn);

%% Simulate the complete identified model with stationary sinusoidal inputs
[sse, y_hat, error, U_delay_mat_sim] = model_simulation_i(model,u_fft,u_fft.*0,20);
figure;plot(y_hat(:,1));

figure;plot( w, abs(Ts.*fft(y_hat(:,1), fftn)) );

%------------------------- Phase/delay analysis ------------------------------------------
%% Decompose the model into canonical low- and high-frequency components
fft_res = Fs/N;
K_LF = round(frq_bndw_LF/fft_res); if K_LF==0; K_LF=1; end %floor( ( ( (frq_bndw_LF)/Fs )*N ) );
K_HF = round(frq_bndw_HF/fft_res); if K_HF==0; K_HF=1; end %floor( ( ( (frq_bndw_HF)/Fs )*N ) );
sig_filt_LF = nrrw_bnd_fft_filt(s_final_trim', Fs, pos_freq_comp(1), K_LF, filt_typ{1});
sig_filt_HF = nrrw_bnd_fft_filt(s_final_trim', Fs, pos_freq_comp(2), K_HF, filt_typ{2});

std_cos = sqrt(0.5);
sig_filt_std = std([sig_filt_LF', sig_filt_HF'] ,1); scl_fctr = sig_filt_std./std_cos;

[LF_sig, AM_sig] = model_simulation_OSD(model,pos_freq_comp,Ts,phi,scl_fctr);
NARX_mod_outpt = LF_sig + AM_sig;

%% Estimate preferred phase from the noise-free canonical model output
LF_phs_mod = angle(hilbert(LF_sig));
HF_env_mod = abs(hilbert(AM_sig));
nbins = 50;
[bin_centers, amp_means] = pac_hist(nbins , LF_phs_mod , HF_env_mod , false);
figure;
bar(rad2deg(bin_centers), amp_means./sum(amp_means), 'FaceColor', 'k', 'EdgeAlpha', 0);
ylabel('Mean Amp.');
set(gca,"TickLabelInterpreter",'latex'); set(gca, 'XTick', [-180 -90 0 90 180]); set(gca, 'XTickLabel', {'$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'});
set(gca,'YTick',[]);
set(gca, 'FontSize', 50);

%% Compare with the known synthetic slow and amplitude-modulated components
LF_phs_mod = angle(hilbert(s_LF));
HF_env_mod = abs(hilbert(s_HF1_shft));
[bin_centers, amp_means] = pac_hist(nbins , LF_phs_mod , HF_env_mod , false);
figure;
bar(rad2deg(bin_centers), amp_means./sum(amp_means), 'FaceColor', [0.47,0.67,0.19], 'EdgeAlpha', 0);
ylabel('Mean Amp.');
set(gca,"TickLabelInterpreter",'latex'); set(gca, 'XTick', [-180 -90 0 90 180]); set(gca, 'XTickLabel', {'$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'});
set(gca,'YTick',[]);
set(gca, 'FontSize', 50);

%% Plot canonical and reference oscillatory components
str = '#ff0000b2'; color_lf = [sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 , 0.7]; % low frequency signal color
str = '#0000ffb2'; color_hf = [sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 , 0.7]; % high freqeuncy singal color
figure; ax1=subplot(1,2,1);plot(LF_sig(1:length(trim_ind)), 'Color', color_lf,'LineWidth',2);axis([0,length(trim_ind),-inf,inf]);
ax2=subplot(1,2,2);plot(AM_sig(1:length(trim_ind)), 'Color', color_hf,'LineWidth',2);axis([0,length(trim_ind),-inf,inf]);
box(ax1,'off');set(get(ax1, 'YAxis'), 'Visible', 'off');set(get(ax1, 'XAxis'), 'Visible', 'off');
box(ax2,'off');set(get(ax2, 'YAxis'), 'Visible', 'off');set(get(ax2, 'XAxis'), 'Visible', 'off');


figure;ax=subplot(1,1,1);plot(LF_sig(1:length(trim_ind)), 'Color', color_lf,'LineWidth',3);axis([0,length(trim_ind),-inf,inf]);hold on;
plot(AM_sig(1:length(trim_ind)), 'Color', color_hf,'LineWidth',2);axis([0,length(trim_ind),-inf,inf]);
box(ax,'off');set(get(ax, 'YAxis'), 'Visible', 'off');set(get(ax, 'XAxis'), 'Visible', 'off'); axis([0,200,-inf,inf]);

figure;ax=subplot(1,1,1);plot(s_LF(trim_ind), 'Color', color_lf,'LineWidth',3);axis([0,length(trim_ind),-inf,inf]);hold on;
plot(s_HF1_shft(trim_ind), 'Color', color_hf,'LineWidth',2);axis([0,length(trim_ind),-inf,inf]);
box(ax,'off');set(get(ax, 'YAxis'), 'Visible', 'off');set(get(ax, 'XAxis'), 'Visible', 'off'); axis([0,200,-inf,inf]);


%% Inspect the canonical spectrum and modulation-index component ratios
MPO_QPC_fft = Ts.*fft(NARX_mod_outpt, fftn);
MPO_QPC_FFT_Mag = abs(MPO_QPC_fft);
MPO_QPC_FFT_Phs = angle(MPO_QPC_fft).*(180/pi);

intrmd_sm = pos_freq_comp(2)+pos_freq_comp(1); intrmd_dff = pos_freq_comp(2)-pos_freq_comp(1);
low_freq_mag  = MPO_QPC_FFT_Mag(freq_pos(pos_freq_comp(1),Fs,25000),:)
high_freq_mag = MPO_QPC_FFT_Mag(freq_pos(pos_freq_comp(2),Fs,25000),:)
intrmd_sm_freq_mag  = MPO_QPC_FFT_Mag(freq_pos(intrmd_sm,Fs,25000),:);
intrmd_dff_freq_mag  = MPO_QPC_FFT_Mag(freq_pos(intrmd_dff,Fs,25000),:);
intrmd_high_freq_ratio = [intrmd_sm_freq_mag/high_freq_mag, intrmd_dff_freq_mag/high_freq_mag]
intrmd_low_freq_ratio = [intrmd_sm_freq_mag/low_freq_mag , intrmd_dff_freq_mag/low_freq_mag]
high_low_freq_ratio = high_freq_mag / low_freq_mag

figure;ax=subplot(2,1,1);plot(w, MPO_QPC_FFT_Mag,  'k', 'LineWidth', 1.5);axis([0,Fs/2,0,inf]);
set(ax,'YTick',[]); box(ax,'off');
subplot(2,1,2);plot(w, MPO_QPC_FFT_Phs);axis([0,Fs/2,0,inf]);

%% Compare time-frequency representations of model and observed signals
tm_frq_plt(NARX_mod_outpt, Fs, fftn);
tm_frq_plt(s_final, Fs, length(s_final));

figure;ax1 = subplot(2,1,1);plot(s_final(trim_ind) , 'Color',  [0.47,0.67,0.19], 'LineWidth', 3); axis([0,200,-inf,inf]);
ax2 = subplot(2,1,2);plot( NARX_mod_outpt(1:length(trim_ind)) , 'Color', 'k', 'LineWidth' , 3); axis([0,200,-inf,inf]);
linkaxes([ax1,ax2], 'x');


%% Test

%% Legacy cluster-extraction consistency check
match_ind_intmod = intr_modul_comp(model);
match_ind_arx = lin_comp(model);
match_ind_QPC_narx = match_ind_intmod | match_ind_arx;

if sum(match_ind_intmod)~=0 && sum(match_ind_arx)~=0
    [~,Y_est_sub] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_QPC_narx);
    tspan_nofrf_trim = tspan_nofrf( length(u_fft) - length(Y_est_sub) + 1:end );
    [~,zc_tspan] = min(abs(tspan_nofrf_trim));

    MPO_QPC_fft = Ts.*fft(Y_est_sub(zc_tspan:end), fftn);
    MPO_QPC_FFT_Mag = abs(MPO_QPC_fft);
    MPO_QPC_FFT_Phs = angle(MPO_QPC_fft).*(180/pi);
else
    MPO_QPC_FFT_Mag = zeros(fftn,1);
    MPO_QPC_FFT_Phs = zeros(fftn,1);
end
%% Local functions
function [bin_centers, amp_means] = pac_hist(nbins,LF_phs,HF_env, plt)
%PAC_HIST Compute mean high-frequency amplitude within slow-phase bins.
%   Inputs are the number of bins, slow instantaneous phase in radians,
%   high-frequency amplitude envelope, and an optional plotting flag.
%   Outputs are phase-bin centres and the corresponding mean amplitudes.
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

%% Test functions
function match_ind = intr_modul_comp(model)
%INTR_MODUL_COMP Select cross-input quadratic terms Sigma_u1u2.
%   Input: MODEL is the selected NonSysID-i model. Output MATCH_IND selects
%   cross-input quadratic terms.
mod_term_char = model{1,16}.Properties.RowNames; % Model term character strings
n_terms = length(mod_term_char); % No. of model terms
pattern = '^u(\d+)\(t-(\d+)\)u(\d+)\(t-(\d+)\)$'; % Define the regular expression pattern
match_ind = zeros(n_terms,1);
for i = 1:n_terms
    term_str = mod_term_char(i); term_str = term_str{1}; % Term identification string
    matches = regexp(term_str, pattern, 'tokens'); % Use regexp to find matches
    % Check if matches are found and omit ux(t-c)ux(t-d)
    if ~isempty(matches)
        % Extract the integer values from the matched tokens
        I_1 = str2double(matches{1}{1});
        I_3 = str2double(matches{1}{3});
        % Check the condition I_1 ? I_3
        if I_1 ~= I_3
            match_ind(i) = 1;
        end
    end
end
if model{11} ~= 0 % Remove bias term index if present
    match_ind = match_ind(1:end-1);
end
end

function match_ind = lin_comp(model)
%LIN_COMP Select all first-order delayed-input terms.
%   Input: MODEL is the selected NonSysID-i model. Output MATCH_IND selects
%   all first-order delayed-input terms.
mod_term_char = model{1,16}.Properties.RowNames; % Model term character strings
n_terms = length(mod_term_char); % No. of model terms
pattern = '^u(\d+)\(t-(\d+)\)$'; % Define the regular expression pattern
match_ind = zeros(n_terms,1);
for i = 1:n_terms
    term_str = mod_term_char(i); term_str = term_str{1}; % Term identification string
    matches = regexp(term_str, pattern, 'match'); % Use regexp to find matches
    % Check if matches are found and omit ux(t-c)ux(t-d)
    if ~isempty(matches)
        match_ind(i) = 1;
    end
end
if model{11} ~= 0 % Remove bias term index if present
    match_ind = match_ind(1:end-1);
end
end

function match_ind = lin_comp_u2(model)
%LIN_COMP_U2 Select first-order high-frequency-input terms Sigma_u2.
%   Input: MODEL is the selected NonSysID-i model. Output MATCH_IND selects
%   high-frequency linear terms.
mod_term_char = model{1,16}.Properties.RowNames; % Model term character strings
n_terms = length(mod_term_char); % No. of model terms
pattern = '^u2\(t-(\d+)\)$'; % Define the regular expression pattern
match_ind = zeros(n_terms,1);
for i = 1:n_terms
    term_str = mod_term_char(i); term_str = term_str{1}; % Term identification string
    matches = regexp(term_str, pattern, 'match'); % Use regexp to find matches
    % Check if matches are found and omit ux(t-c)ux(t-d)
    if ~isempty(matches)
        match_ind(i) = 1;
    end
end
if model{11} ~= 0 % Remove bias term index if present
    match_ind = match_ind(1:end-1);
end
end

function match_ind = lin_comp_u1(model)
%LIN_COMP_U1 Select first-order low-frequency-input terms Sigma_u1.
%   Input: MODEL is the selected NonSysID-i model. Output MATCH_IND selects
%   low-frequency linear terms.
mod_term_char = model{1,16}.Properties.RowNames; % Model term character strings
n_terms = length(mod_term_char); % No. of model terms
pattern = '^u1\(t-(\d+)\)$'; % Define the regular expression pattern
match_ind = zeros(n_terms,1);
for i = 1:n_terms
    term_str = mod_term_char(i); term_str = term_str{1}; % Term identification string
    matches = regexp(term_str, pattern, 'match'); % Use regexp to find matches
    % Check if matches are found and omit ux(t-c)ux(t-d)
    if ~isempty(matches)
        match_ind(i) = 1;
    end
end
if model{11} ~= 0 % Remove bias term index if present
    match_ind = match_ind(1:end-1);
end
end
