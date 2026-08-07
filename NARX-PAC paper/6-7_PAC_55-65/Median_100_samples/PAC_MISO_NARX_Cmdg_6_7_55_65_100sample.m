%PAC_MISO_NARX_CMDG_6_7_55_65_100SAMPLE Generate the NARX-PAC ensemble used for the Figure 23 repeatability study.
%   One hundred pink-noise realisations are used to form a 6-7 Hz slow
%   oscillation coupled to a 55-60 Hz fast oscillation. The script stores
%   the raw comodulograms, discriminator maps, model data, and frequency grids.
%
clear;clc;close all;

addpath('\<path-to>\NonSysID-i\');
addpath('\<path-to>\NARX_PAC\');
addpath('\<path-to>\NARX_PAC\Utils\');
%% Set sampling and plotting parameters
Fs = 1000; Ts = 1/Fs;
R=4;C=1;
%% Define numerical helper functions
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);

%% =====================================================
%% PAC LF-sine HF-sine simple model
%% =====================================================

%% Configure the analysis interval
tspan = 0:Ts:25;%(N*Ts-Ts);
N = length(tspan);
fftn = 4000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);

%% Sys ID params

filt_typ = {'sbp','sbp'}; % bw , sbp , guss
frq_bndw_LF = 1;
frq_bndw_HF = 0.5;
disp(['frq_bndw_LF = ', num2str(frq_bndw_LF), ', frq_bndw_HF = ', num2str(frq_bndw_HF)]);

fL_vals = [2:1:15];
fH_vals = [20:1:100];

%Down-sample params
dwn_smpl_F = 250;
dwn_smpl_T = 1/dwn_smpl_F;
dwn_smpl = Fs/dwn_smpl_F;

%% PAC signal params

SNR = 3;
tm_windw = 10;

n_smpls = 100;
rng(200,"twister"); rng_seeds = randi([1,1e5],n_smpls,1);

save_data = false;

% PAC model params
LF_freq_1 = [6,7]; HF_freq_1 = [55,60];
disp(['LF = ', num2str(LF_freq_1), ', HF1 = ', num2str(HF_freq_1)]);
m = 0.5; A = 200; C = 1*1e-6; am_lag = 16;

% PAC signal trimming params
tm_itr = 0;
trim_ind = (tm_itr*tm_windw/dwn_smpl_T)+1:((tm_itr+1)*tm_windw)/dwn_smpl_T;%(1 + (30/fL)/Ts);%
len_trim_ind = length(trim_ind);

%% Nonsine PAC single LF coupling a single HF

s_final_mat = zeros(len_trim_ind, n_smpls);

for i = 1:n_smpls
    Fs = 1000; Ts = 1/Fs;
    tspan = 0:Ts:25;%(N*Ts-Ts);
    N = length(tspan);

    %=======================
    % Low- and high-frequency oscillations from random pink noise.

    rng(rng_seeds(i), 'twister'); [s_LF, s_HF, ~, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_1);
    rng(rng_seeds(i)+1000, 'twister'); [~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, LF_freq_1, HF_freq_1);
    %=======================
    % Complex PAC

    [s_final, s_HF1_shft] = pac_general_1(s_LF, s_HF, A, C, m, am_lag);

    if am_lag~=0
        s_final = s_final(am_lag:end);
        tspan = tspan(am_lag:end);
        N = length(tspan);
    end
    %=======================
    %Down sample

    s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
    s_final = s_final(1:dwn_smpl:N);
    tspan = tspan(1:dwn_smpl:N);
    Fs = dwn_smpl_F; Ts = 1/Fs;

    pink = pink(1:dwn_smpl:N);
    %=======================

    %=======================
    %Trim the PAC signal to get a small segment

    s_final_trim = s_final(trim_ind)';
    N = length(s_final_trim);
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
    %% Complete the ensemble loop

end

%% NARX based MISO PAC Comodulogram
Comods = cell(1,n_smpls);
Comods_diff = cell(1,n_smpls);
All_freq_comb_dat = cell(1,n_smpls);
PAC_data_mat = cell(1,n_smpls);

for j = 1:n_smpls
    %=======================
    [Comod , diff_comod, phs_data_mat, fL_grd, fH_grd, All_freq_comb_1, All_freq_comb, All_freq_comb_ARX_1, All_freq_comb_ARX_2, narx_pac_modls_1 , narx_pac_modls_2]...
        = pac_miso_Cmdg_mod_21(s_final_mat(:,j), fL_vals, fH_vals, Fs, 3, filt_typ, frq_bndw_LF, frq_bndw_HF);

    Comods{1,j}            = Comod;
    Comods_diff{1,j}       = diff_comod;
    All_freq_comb_dat{1,j} = All_freq_comb_1;
    PAC_data_mat{1,j}      = phs_data_mat;

    disp(j);

    %=======================
end

%% Save the ensemble results
if save_data
    file_dir = '/home/gunawardes/Documents/Matlab/CFC/Results/MISO_NARX/SyntheticData/Multisample/1/';
    file_name = ['F2_', num2str(SNR) ,'_', num2str(tm_windw) ,'s_wrk.mat'];
    save([file_dir,file_name]);
end
%% Run post-processing after this
% Run the file PAC_MISO_NARX_mltsmpl_dat_prcss.m
%% Local functions - PAC general

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



