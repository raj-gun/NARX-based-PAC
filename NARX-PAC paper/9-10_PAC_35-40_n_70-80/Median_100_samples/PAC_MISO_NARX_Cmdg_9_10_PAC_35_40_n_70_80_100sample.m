clear all;clc;close all;

addpath('\<path-to>\NonSysID-i\');
addpath('\<path-to>\NARX_PAC\');
addpath('\<path-to>\NARX_PAC\Utils\');
%%
Fs = 1000; Ts = 1/Fs;
R=4;C=1;
%%
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);

%% =====================================================
%% PAC LF-sine HF-sine simple model
%% =====================================================

%%
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
rng(100,"twister"); rng_seeds = randi([1,1e5],n_smpls,1);

save_data = false;

LF_freq_1 = [9,10];
HF_freq_1 = [65,75]+5; HF_freq_2 = [40,45]-5; % Working

disp(['LF = ', num2str(LF_freq_1), ', HF1 = ', num2str(HF_freq_1), ', HF2 = ', num2str(HF_freq_2)]);

m1 = 0.5; A1 = 50;
m2 = 0.25; A2 = 150;
am_lag_1 = 16; am_lag_2 = 90; am_lag = max([am_lag_1,am_lag_2]);

%PAC signal trimming params
tm_itr = 0;
trim_ind = (tm_itr*tm_windw/dwn_smpl_T)+1:((tm_itr+1)*tm_windw)/dwn_smpl_T;%(1 + (30/fL)/Ts);%
len_trim_ind = length(trim_ind);

%% Nonsine PAC single LF coupling with two distinct HFs

s_final_mat = zeros(len_trim_ind, n_smpls);

for i = 1:n_smpls
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

    s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
    s_final = s_final(1:dwn_smpl:N);
    tspan   = tspan(1:dwn_smpl:N);
    Fs      = dwn_smpl_F; Ts = 1/Fs;

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

    fL_diff      = mean(abs(diff(fL_vals))); fH_diff = mean(abs(diff(fH_vals)));
    rect_pos_box = @(LF_freq, HF_freq, fL_diff, fH_diff) [LF_freq(1)-fL_diff*0.5, HF_freq(1)-fH_diff*0.5, (abs(diff(LF_freq))*fL_diff)+1, (abs(diff(HF_freq))*1)+1];
    rect_pos_1   = rect_pos_box(LF_freq_1, HF_freq_1, fL_diff, fH_diff);
    rect_pos_2   = rect_pos_box(LF_freq_1, HF_freq_2, fL_diff, fH_diff);

    disp(j);

    %=======================
end

%%
if save_data
    file_dir = '/home/gunawardes/Documents/Matlab/CFC/Results/MISO_NARX/SyntheticData/Multisample/1/';
    file_name = ['F1_', num2str(SNR) ,'_', num2str(tm_windw) ,'s_wrk.mat'];
    save([file_dir,file_name]);
end
%%
%% Run post-processing after this
% Run the file PAC_MISO_NARX_mltsmpl_dat_prcss.m
%%

%% Local functions - PAC general

% Equation adapted from Jiang et al., (2015) NeuroImage
function [s_final, s_HF1_shft] = pac_general_1(s_LF, s_HF, a, c, m, delay_ind)
s_HF1 = m .* ( 1 - ( 1./(1 + exp(-a.*(s_LF-c))) ) ) .* s_HF;
if delay_ind ~= 0; s_HF1_shft = zeros(size(s_HF1)); s_HF1_shft(delay_ind:end) = s_HF1(1:end-delay_ind+1);
else; s_HF1_shft = s_HF1; end
s_final = s_LF + s_HF1_shft;
end
%-----------------------------------------------------

% Simplest form of PAC generaltion in electronics
%(J. Smith, Mathematics of the discrete Fourier transform (DFT). [North Charleston]: BookSurge, 2010.)
function [s_final, s_HF1_shft] = pac_simple(s_LF, s_HF, a, m, delay_ind)
s_HF1 = m .* (1 + a.*s_LF) .* s_HF;
if delay_ind~=0; s_HF1_shft = zeros(size(s_HF1)); s_HF1_shft(delay_ind:end) = s_HF1(1:end-delay_ind+1); else; s_HF1_shft = s_HF1; end
s_final = s_LF + s_HF1_shft;
%%
end



