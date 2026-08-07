clear all;clc;close all;
%% Load data

file_dir = '\<path-to>\NARX-PAC paper\LFP\Plot data\';
file_name = '500Hz_200s-0wndw_1-0.5-sbp_sbp_LFcos_2.mat';
MISO_NARX_HG = load([file_dir,file_name]);

file_name = '500Hz_200s-0wndw_1-0.5-sbp_sbp_LFcos.mat';
MISO_NARX_HFO = load([file_dir,file_name]);


file_dir = '\<path-to>\NARX-PAC paper\LFP\Plot data\';
file_name = 'HFO.mat';
Othr_Mthds_HFO = load([file_dir,file_name]);

file_dir = '\<path-to>\NARX-PAC paper\LFP\Plot data\';
file_name = 'HG.mat';
Othr_Mthds_HG = load([file_dir,file_name]);
% plot_data = { {OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins} , ...
%               {GLM_org, GLM_2, GLM_robust, flow_GLM, fhigh_GLM} };

%% Plots PAC


figure;
tiledlayout(2,4, "TileSpacing", 'loose');
font_size = 16;

% Tile 1
nexttile
imagesc(Othr_Mthds_HFO.plot_data{1,1}{1,5}, Othr_Mthds_HFO.plot_data{1,1}{1,6}, Othr_Mthds_HFO.plot_data{1,1}{1,1}');  colorbar; axis xy; set(gca, 'FontSize', font_size);
title({'Ozkurt et. al. 2011';' '});
ylabel({'\bf HFO' ; ' ' ; '\rm High Frequency (Hz)'});
% Tile 2
nexttile
imagesc(Othr_Mthds_HFO.plot_data{1,1}{1,5}, Othr_Mthds_HFO.plot_data{1,1}{1,6}, Othr_Mthds_HFO.plot_data{1,1}{1,2}');  colorbar; axis xy; set(gca, 'FontSize', font_size);
title({'Canolty et. al. 2010';' '});
% Tile 3
nexttile
imagesc(Othr_Mthds_HFO.plot_data{1,2}{1,4}, Othr_Mthds_HFO.plot_data{1,2}{1,5}, Othr_Mthds_HFO.plot_data{1,2}{1,3}');  colorbar; axis xy; set(gca, 'FontSize', font_size);
title({'Penny et. al. 2008';' '});
% Tile 4
nexttile
imagesc(MISO_NARX_HFO.fL_vals, MISO_NARX_HFO.fH_vals, MISO_NARX_HFO.Comods{1});  colorbar; axis xy; set(gca, 'FontSize', font_size);
title({'NARX-based PAC';' '});

% Tile 5
nexttile
imagesc(Othr_Mthds_HG.plot_data{1,1}{1,5}, Othr_Mthds_HG.plot_data{1,1}{1,6}, Othr_Mthds_HG.plot_data{1,1}{1,1}');  colorbar; axis xy; set(gca, 'FontSize', font_size);
xlabel('Low Frequency (Hz)');
ylabel({'\bf HG' ; ' ' ; '\rm High Frequency (Hz)'});
% Tile 6
nexttile
imagesc(Othr_Mthds_HG.plot_data{1,1}{1,5}, Othr_Mthds_HG.plot_data{1,1}{1,6}, Othr_Mthds_HG.plot_data{1,1}{1,2}');  colorbar; axis xy; set(gca, 'FontSize', font_size);
xlabel('Low Frequency (Hz)');
% Tile 7
nexttile
imagesc(Othr_Mthds_HG.plot_data{1,2}{1,4}, Othr_Mthds_HG.plot_data{1,2}{1,5}, Othr_Mthds_HG.plot_data{1,2}{1,3}');  colorbar; axis xy; set(gca, 'FontSize', font_size);
xlabel('Low Frequency (Hz)');
% Tile 8
nexttile
imagesc(MISO_NARX_HG.fL_vals, MISO_NARX_HG.fH_vals, MISO_NARX_HG.Comods{1});  colorbar; axis xy; set(gca, 'FontSize', font_size);
xlabel('Low Frequency (Hz)');

%% PLots Phs data
file_dir = '\<path-to>\NARX-PAC paper\LFP\Plot data\';
file_name = 'phs_PltDat.mat';
MISO_NARX_phs_HFO = load([file_dir,file_name]);

file_dir = '\<path-to>\NARX-PAC paper\LFP\Plot data\';
file_name = 'phs_PltDat_2.mat';
MISO_NARX_phs_HG = load([file_dir,file_name]);

MISO_NARX_phs_HFO_ind = find(MISO_NARX_phs_HFO.plot_data{1,2}==8);
MISO_NARX_phs_HG_ind = find(MISO_NARX_phs_HG.plot_data{1,2}==8);

HF_freqs_HFO = MISO_NARX_phs_HFO.plot_data{1,1}{MISO_NARX_phs_HFO_ind}{1,2};
HF_freqs_HG = MISO_NARX_phs_HG.plot_data{1,1}{MISO_NARX_phs_HG_ind}{1,2};

high_freq_HFO = 33;
high_freq_HG = 80;
%% PLots Phs
%% PLots Phs 2
figure('Renderer','Painters');
tiledlayout(2,4, 'TileSpacing', 'compact', 'Padding','compact');
font_size = 12;

% Tile 1
nexttile
surf(MISO_NARX_phs_HFO.plot_data{1,1}{MISO_NARX_phs_HFO_ind}{1,1}' , MISO_NARX_phs_HFO.plot_data{1,1}{MISO_NARX_phs_HFO_ind}{1,5} , MISO_NARX_phs_HFO.plot_data{1,1}{MISO_NARX_phs_HFO_ind}{1,3}','EdgeColor','none'); 
set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
axis([-pi,pi, min(HF_freqs_HFO),max(HF_freqs_HFO), -inf,inf]);
set(gca, 'FontSize', font_size); view(2); colorbar;
ylabel({'\bf HFO' ; ' ' ; '\rm High Frequency (Hz)'});
title({'NARX-based PAC';' '});

% Tile 2
nexttile
Phs_Amp = Othr_Mthds_HFO.plot_data{1,1}{1,4};
phs_bins = Othr_Mthds_HFO.plot_data{1,1}{1,7};
flow_MI = Othr_Mthds_HFO.plot_data{1,1}{1,5};
fhigh_MI = Othr_Mthds_HFO.plot_data{1,1}{1,6};
phs_freq = reshape(Phs_Amp(flow_MI==8,:,:), size(Phs_Amp,2), size(Phs_Amp,3) );
[~,idxsIntoA] = intersect(fhigh_MI,HF_freqs_HFO,'stable');
fhigh_MI = fhigh_MI(idxsIntoA);
phs_freq = phs_freq(idxsIntoA,:);
imagesc(phs_bins, fhigh_MI, phs_freq ); axis xy;
set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
axis([-pi,pi, min(HF_freqs_HFO),max(HF_freqs_HFO), -inf,inf]);
set(gca, 'FontSize', font_size); colorbar;
title({'Filtering-based';' '});

%Tile 1
nexttile
HF_freq_ind = find(HF_freqs_HFO == high_freq_HFO);
bar(MISO_NARX_phs_HFO.plot_data{1,1}{MISO_NARX_phs_HFO_ind}{1,1}(:,HF_freq_ind), MISO_NARX_phs_HFO.plot_data{1,1}{MISO_NARX_phs_HFO_ind}{1,3}(:,HF_freq_ind),  'FaceColor', [0.2 0.6 0.8], 'EdgeAlpha', 0);
set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
set(gca, 'FontSize', font_size);
ylabel({'\bf HFO' ; ' ' ; ['\rm Mean Ampl. at ', num2str(high_freq_HFO), 'Hz']}); %ylabel(['\rm Mean Ampl. at ', num2str(high_freq_HFO), 'Hz']);
title({'NARX-based PAC';' '});

%Tile 2
nexttile
Phs_Amp = Othr_Mthds_HFO.plot_data{1,1}{1,4};
phs_bins = Othr_Mthds_HFO.plot_data{1,1}{1,7};
flow_MI = Othr_Mthds_HFO.plot_data{1,1}{1,5};
fhigh_MI = Othr_Mthds_HFO.plot_data{1,1}{1,6};
phs_freq = reshape(Phs_Amp(flow_MI==8,:,:), size(Phs_Amp,2), size(Phs_Amp,3) );
[~,idxsIntoA] = intersect(fhigh_MI,HF_freqs_HFO,'stable');
fhigh_MI = fhigh_MI(idxsIntoA);
phs_freq = phs_freq(idxsIntoA,:);
HF_freq_ind = find(fhigh_MI == high_freq_HFO);
bar(phs_bins, phs_freq(HF_freq_ind,:),  'FaceColor', [0.2 0.6 0.8], 'EdgeAlpha', 0);
set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
set(gca, 'FontSize', font_size);
title({'Filtering-based';' '});

% Tile 3
nexttile
surf(MISO_NARX_phs_HG.plot_data{1,1}{MISO_NARX_phs_HG_ind}{1,1}' , MISO_NARX_phs_HG.plot_data{1,1}{MISO_NARX_phs_HG_ind}{1,5} , MISO_NARX_phs_HG.plot_data{1,1}{MISO_NARX_phs_HG_ind}{1,3}','EdgeColor','none'); 
set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
axis([-pi,pi, min(HF_freqs_HG),max(HF_freqs_HG), -inf,inf]);
set(gca, 'FontSize', font_size); view(2); colorbar;
xlabel('Low Frequency Phase (rads)');
ylabel({'\bf HG' ; ' ' ; '\rm High Frequency (Hz)'});

% Tile 4
nexttile
Phs_Amp = Othr_Mthds_HG.plot_data{1,1}{1,4};
phs_bins = Othr_Mthds_HG.plot_data{1,1}{1,7};
flow_MI = Othr_Mthds_HG.plot_data{1,1}{1,5};
fhigh_MI = Othr_Mthds_HG.plot_data{1,1}{1,6};
phs_freq = reshape(Phs_Amp(flow_MI==8,:,:), size(Phs_Amp,2), size(Phs_Amp,3) );
[~,idxsIntoA] = intersect(fhigh_MI,HF_freqs_HG,'stable');
fhigh_MI = fhigh_MI(idxsIntoA);
phs_freq = phs_freq(idxsIntoA,:);
imagesc(phs_bins, fhigh_MI, phs_freq ); axis xy;
set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
axis([-pi,pi, min(HF_freqs_HG),max(HF_freqs_HG), -inf,inf]); 
set(gca, 'FontSize', font_size); colorbar;
xlabel('Low Frequency Phase (rads)');

%Tile 3
nexttile
HF_freq_ind = find(HF_freqs_HG == high_freq_HG);
bar(MISO_NARX_phs_HG.plot_data{1,1}{MISO_NARX_phs_HG_ind}{1,1}(:,HF_freq_ind), MISO_NARX_phs_HG.plot_data{1,1}{MISO_NARX_phs_HG_ind}{1,3}(:,HF_freq_ind),  'FaceColor', [0.2 0.6 0.8], 'EdgeAlpha', 0);
set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
set(gca, 'FontSize', font_size);
xlabel('Low Frequency Phase (rads)');
ylabel({'\bf HG' ; ' ' ; ['\rm Mean Ampl. at ', num2str(high_freq_HG), 'Hz']}); %ylabel(['\rm Mean Ampl. at ', num2str(high_freq_HG), 'Hz']);

%Tile 4
nexttile
Phs_Amp = Othr_Mthds_HG.plot_data{1,1}{1,4};
phs_bins = Othr_Mthds_HG.plot_data{1,1}{1,7};
flow_MI = Othr_Mthds_HG.plot_data{1,1}{1,5};
fhigh_MI = Othr_Mthds_HG.plot_data{1,1}{1,6};
phs_freq = reshape(Phs_Amp(flow_MI==8,:,:), size(Phs_Amp,2), size(Phs_Amp,3) );
[~,idxsIntoA] = intersect(fhigh_MI,HF_freqs_HG,'stable');
fhigh_MI = fhigh_MI(idxsIntoA);
phs_freq = phs_freq(idxsIntoA,:);
HF_freq_ind = find(fhigh_MI == high_freq_HG);
bar(phs_bins, phs_freq(HF_freq_ind,:),  'FaceColor', [0.2 0.6 0.8], 'EdgeAlpha', 0);
set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
set(gca, 'FontSize', font_size);
xlabel('Low Frequency Phase (rads)');
%% PLots Phs 3







