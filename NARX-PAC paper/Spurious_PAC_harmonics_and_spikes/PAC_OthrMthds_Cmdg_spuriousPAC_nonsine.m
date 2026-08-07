%PAC_OTHRMTHDS_CMDG_SPURIOUSPAC_NONSINE Analyse a non-sinusoidal waveform using conventional PAC methods.
%   The script calculates the Ozkurt, Canolty, Tort, and GLM results used in
%   Figure 20 to compare harmonic-related false detections with the
%   post-processed NARX-PAC result.
%
clear all;clc;close all;

addpath('\<path-to>\NARX_PAC\Utils\');
addpath('\<path-to>\Other_Methods\');
%% Set sampling and plotting parameters
Fs = 1000; Ts = 1/Fs;
R=4;C=1;
%% Define numerical helper functions
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);

%% =====================================================
%% Configure the spurious-PAC test signal
%% =====================================================

%% Configure the analysis interval
tspan = 0:Ts:25-Ts;%(N*Ts-Ts);
N = length(tspan);
fftn = 4000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);

%% Generate the non-sinusoidal waveform

LF = 7;
u = sin( (2*pi*LF).*tspan )./LF;
s_final = 10 ./ ( 1 + exp( -12.*( 0.5.*(1+5.*u) -0.7 ) )  );
am_lag = 0;


n_smpls = 100; rng(100,"twister"); rng_seeds = randi([1,1e5],n_smpls,1); 
rng( rng_seeds(60) ); [~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, [9.5,10.5], [55,65]);

%% Downsample the signal

% 


%% Select and trim one signal segment

%Trim the PAC signal to get a small segment
tm_windw = 20; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
s_final_trim = s_final(trim_ind)';%(600:1100)';

rng(100,"twister");
SNR = 3;
pink = pink(trim_ind)';
s_final_trim = ( SNR*( std(pink)/std(s_final_trim) ) ) .* s_final_trim;
sn_ratio = snr(s_final_trim,pink); disp(['SNR = ', num2str(sn_ratio), 'dB or ',  num2str(db2mag(sn_ratio))]);
s_final_trim = pink + s_final_trim;

N = length(s_final_trim);
figure;plot(tspan(trim_ind),s_final_trim);
%----------------------
fftn = length(s_final_trim); 
w = 0:Fs/fftn:Fs-(Fs/fftn);
tm_frq_plt(s_final_trim, Fs, fftn);

fL_vals = [1:1:20]; fH_vals = [15:1:90];


%% Compute conventional PAC comodulograms
[OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins] = modulationindex_directestimate_mod(s_final_trim', Fs,fL_vals,fH_vals,0.5,0.5,100);

phs_freq = reshape(Phs_Amp(flow_MI==6,:,:), size(Phs_Amp,2), size(Phs_Amp,3) );
figure;imagesc(phs_bins, fhigh_MI, phs_freq ); axis xy;

figure
imagesc(flow_MI, fhigh_MI, OzktMI');  colorbar;
axis xy
set(gca, 'FontSize', 18);

figure
imagesc(flow_MI, fhigh_MI, CanltyMI');  colorbar;
axis xy
set(gca, 'FontSize', 18);

figure
imagesc(flow_MI, fhigh_MI, TortMI');  colorbar;
axis xy
set(gca, 'FontSize', 18);

[GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM] = general_linear_index_mod(s_final_trim', Fs,fL_vals,fH_vals,0.5,0.5);

figure
imagesc(flow_GLM, fhigh_GLM, GLM_org');  colorbar;
axis xy
set(gca, 'FontSize', 18);


figure;
imagesc(flow_GLM, fhigh_GLM, GLM_2');  colorbar;
axis xy
set(gca, 'FontSize', 18);


figure
imagesc(flow_GLM, fhigh_GLM, GLM_robust');  colorbar;
axis xy
set(gca, 'FontSize', 18);

plot_data = { {OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins} , ...
              {GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM} };

% ===============================
% Use the following code snippet to save the workspace for plotting

% save_file_name = 'NonSine';
% save_file_dir = '\<path-to>\NARX-PAC paper\Spurious_PAC_harmonics_and_spikes\Data_Other_Methods\';
% save([save_file_dir, save_file_name], 'plot_data');
% ===============================

%% =====================================================
%% Local functions
%% =====================================================

%% Van der Pol waveform
function [s_LF] = van_d_pol_LF(tspan,w0)
%VAN_D_POL_LF Generate a normalised non-sinusoidal Van der Pol waveform.
%   TSPAN is the integration time vector, and W0 sets the oscillator angular
%   frequency. S_LF is the resulting zero-mean, unit-peak waveform.
ep=5;%w0=100;

dEqs = @(t, x) [
    x(2);
    ep*w0*(1-x(1)^2)*x(2) - x(1)*w0^2;
    ];
% Solve the system of differential equations
[~,x] = ode45(dEqs, tspan, [2 1]);
s_LF = x(:,1)./max(x(:,1));
s_LF = s_LF - mean(s_LF);

end