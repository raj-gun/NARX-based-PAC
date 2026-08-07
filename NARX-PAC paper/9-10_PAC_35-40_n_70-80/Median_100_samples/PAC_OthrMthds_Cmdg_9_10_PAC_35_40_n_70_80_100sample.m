%PAC_OTHRMTHDS_CMDG_9_10_PAC_35_40_N_70_80_100SAMPLE Generate conventional PAC results for the comparisons in Figures 24-26.
%   The two-band synthetic PAC realisations are analysed using the Ozkurt,
%   Canolty, Tort, and GLM measures. The resulting arrays remain in the
%   workspace for subsequent median-comodulogram plots.
%
clear;clc;close all;

addpath('\<path-to>\NARX_PAC\Utils\');
addpath('\<path-to>\Methods\Matlab_Code');
%% Set sampling and plotting parameters
Fs = 1000; Ts = 1/Fs;
R=4;C=1;
%% Define numerical helper functions
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);

%% =====================================================
%% Configure the synthetic PAC experiment
%% =====================================================

%% Configure the analysis interval
tspan = 0:Ts:25;%(N*Ts-Ts);
N = length(tspan);
fftn = 4000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);

%% Set system-identification parameters

filt_typ = {'sbp','sbp'}; % bw , sbp , guss
frq_bndw_LF = 1;
frq_bndw_HF = 0.5;
disp(['frq_bndw_LF = ', num2str(frq_bndw_LF), ', frq_bndw_HF = ', num2str(frq_bndw_HF)]);

fL_vals = [2:1:15];
fH_vals = [20:1:100];

% % Downsampling parameters

%% Set PAC signal parameters

SNR = 3;
n_smpls = 100; 
rng(100,"twister"); rng_seeds = randi([1,1e5],n_smpls,1);

LF_freq_1 = [9,10];
HF_freq_1 = [65,75]+5; HF_freq_2 = [40,45]-5; % Working

disp(['LF = ', num2str(LF_freq_1), ', HF1 = ', num2str(HF_freq_1), ', HF2 = ', num2str(HF_freq_2)]);

m1 = 0.5; A1 = 50;
m2 = 0.25; A2 = 150;
am_lag_1 = 16; am_lag_2 = 90; am_lag = max([am_lag_1,am_lag_2]);

% PAC signal trimming parameters
tm_windw = 10; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
len_trim_ind = length(trim_ind);

%% Generate ensemble PAC signals

s_final_mat = zeros(len_trim_ind, n_smpls);

parfor i = 1:n_smpls
    Fs = 1000; Ts = 1/Fs;
    tspan = 0:Ts:25;%(N*Ts-Ts);
    N = length(tspan);

    %=======================
    % Low- and high-frequency oscillations from random pink noise.
    
    rng(rng_seeds(i), 'twister'); [s_LF_1, s_HF_1, pink, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_1);
    [~, s_HF_2, ~, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_2);
    rng(rng_seeds(i)+1000, 'twister'); [~, ~, pink1, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_2);
    %=======================
    % Simple PAC

    [~, s_HF1_shft_1] = pac_simple(s_LF_1, s_HF_1, A1, m1, am_lag_1);
    [~, s_HF1_shft_2] = pac_simple(s_LF_1, s_HF_2, A2, m2, am_lag_2);
    s_final = s_LF_1 +  s_HF1_shft_1 + s_HF1_shft_2;
    %=======================
    if am_lag~=0
        s_final = s_final(am_lag:end);
        tspan = tspan(am_lag:end);
        N = length(tspan);
    end
    %=======================
    %Down sample


    %=======================

    %=======================
    %Trim the PAC signal to get a small segment

    s_final_trim = s_final(trim_ind)';
    %----------------------
    fftn = length(s_final_trim);
    %=======================

    %=======================
    % Add pink noise

    pink         = pink(trim_ind)';
    s_final_trim = ( SNR*( std(pink)/std(s_final_trim) ) ) .* s_final_trim;
    s_final_trim = pink + s_final_trim;
    %=======================

    s_final_mat(:,i) = s_final_trim;

end

%% Compute conventional PAC comodulograms

h = waitbar(0, ['Dataset ',0,'/',num2str(n_smpls)]);
Comods = cell(7,n_smpls);
for i = 1:n_smpls
    
    [OzktMI, CanltyMI, TortMI, Phs_Amp, flow_MI, fhigh_MI, phs_bins] = modulationindex_directestimate_mod(s_final_mat(:,i)', Fs,fL_vals,fH_vals,0.5,0.5,100);
    [~, ~, GLM_robust, flow_GLM, fhigh_GLM] = general_linear_index_mod(s_final_mat(:,i)', Fs,fL_vals,fH_vals,0.5,0.5);
    
    Comods{1,i} = OzktMI;
    Comods{2,i} = CanltyMI;
    Comods{3,i} = TortMI;
    Comods{4,i} = Phs_Amp;
    Comods{5,i} = flow_MI;
    Comods{6,i} = fhigh_MI;
    Comods{7,i} = phs_bins;
    Comods{8,i} = GLM_robust;
    Comods{9,i} = flow_GLM;
    Comods{10,i} = fhigh_GLM;
    
    h = waitbar(i/n_smpls, h, ['Dataset ',num2str(i),'/',num2str(n_smpls)]);
end

Comods_mat_OzktMI = cat(3, Comods{1,:});
Comods_mat_CanltyMVL = cat(3, Comods{2,:});
Comods_mat_TortMI = cat(3, Comods{3,:});
Comods_mat_GLM_robust = cat(3, Comods{8,:});

Comods{5,1} = fL_vals;
Comods{6,1} = fH_vals;
%% Plot the ensemble conventional-method comodulograms


figure; tiledlayout(2,2);
nexttile;
imagesc(fL_vals, fH_vals, mean(Comods_mat_OzktMI,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rect_pos_2 = [LF_freq_1(1)-0.25, HF_freq_2(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_2))*1)+1]; rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);

nexttile;
imagesc(fL_vals, fH_vals, mean(Comods_mat_CanltyMVL,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rect_pos_2 = [LF_freq_1(1)-0.25, HF_freq_2(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_2))*1)+1]; rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);

nexttile;
imagesc(fL_vals, fH_vals, mean(Comods_mat_TortMI,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rect_pos_2 = [LF_freq_1(1)-0.25, HF_freq_2(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_2))*1)+1]; rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);

nexttile;
imagesc(flow_GLM, fhigh_GLM, mean(Comods_mat_GLM_robust,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
rect_pos_2 = [LF_freq_1(1)-0.25, HF_freq_2(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_2))*1)+1]; rectangle('Position',rect_pos_2, 'EdgeColor','r', 'LineWidth', 1);



%% Local PAC signal-generation functions

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



