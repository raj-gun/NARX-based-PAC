%PAC_OTHRMTHDS_CMDG_TORTLABLFP2 Apply conventional PAC methods to the Tortlab high-gamma LFP recording.
%   The 200 s recording segment is analysed using the Ozkurt, Canolty, Tort,
%   and GLM measures. Results are stored in PLOT_DATA for comparison with the
%   NARX-PAC results.
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
%% Load the LFP recording

load('\<path-to>\NARX-PAC paper\LFP\LFP data\LFP_HG_HFO.mat');
s_final = lfpHG;
%or
N = length(s_final);
tspan = 0:Ts:(N*Ts-Ts);
fftn = 4000;%Fs/N;
tm_frq_plt(s_final, Fs, fftn); % Visualise PAC  in time-frequency plots
%% Downsample the LFP signal

% 


%% Select the analysis segment
tm_windw = 200; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
s_final_trim = s_final(trim_ind)';%(600:1100)';
N = length(s_final_trim);
figure;plot(tspan(trim_ind),s_final_trim);
%----------------------
fftn = length(s_final_trim); 
w = 0:Fs/fftn:Fs-(Fs/fftn);
tm_frq_plt(s_final_trim, Fs, fftn);

fL_vals = [3:1:13]; fH_vals = [30:1:200];


%% Compute conventional PAC comodulograms
[OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins] = modulationindex_directestimate_mod(s_final_trim', Fs,fL_vals,fH_vals,0.5,0.5,100);

phs_freq = reshape(Phs_Amp(flow_MI==8,:,:), size(Phs_Amp,2), size(Phs_Amp,3) );
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

% Save plot_data as HG.mat