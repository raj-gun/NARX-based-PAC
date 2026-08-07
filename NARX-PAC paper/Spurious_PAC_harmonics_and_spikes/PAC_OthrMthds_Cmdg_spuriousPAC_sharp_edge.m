%PAC_OTHRMTHDS_CMDG_SPURIOUSPAC_SHARP_EDGE Analyse a sharp-edged waveform using conventional PAC methods.
%   The script calculates the Ozkurt, Canolty, Tort, and GLM results used in
%   Figure 21 to compare harmonic-related false detections with the
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

%% Generate the sharp-edged waveform
LF = 9;
s_final = sharp_edge(LF, 0.2, Ts, N);

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

% save_file_name = 'Sharp_Edge';
% save_file_dir = '\<path-to>\NARX-PAC paper\Spurious_PAC_harmonics_and_spikes\Data_Other_Methods\';
% save([save_file_dir, save_file_name], 'plot_data');
% ===============================

%% =====================================================
%% Local functions
%% =====================================================

%% Sharp-edged waveform
% Code adapted from Kramer et al. (2008) and Ozkurt et al. (2011),
% Journal of Neuroscience Methods.
function [s] = sharp_edge(f, edge_pos, Ts, N)
%SHARP_EDGE Generate a periodic waveform with an abrupt edge.
%   F is the base frequency, EDGE_POS specifies the cut position within each
%   period, TS is the sampling interval, and N is the requested sample count.
%   S is the resulting sharp-edged waveform.

T = 1/f;
cut_point = round(edge_pos*T/Ts);
join_point = cut_point + round(0.1*T/Ts);

len_s1 = T/Ts;
len_s2 = floor(len_s1-(join_point-cut_point-1));
n_cycles = floor(N/len_s2);
s = zeros(1,N);

for k=1:n_cycles
    s1 = 2.*cos(2.*pi.*(0:Ts*f:1-Ts*f));
    s2 = zeros(1,len_s2);
    
    s2(1:cut_point) = s1(1:cut_point);
    s2(cut_point+1:end) = s1(join_point:end);

    s( (k-1)*len_s2 + 1 : k*len_s2  ) = s2;
end

end
