clear;clc;close all;

addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\Utils\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\Methods\Matlab_Code');
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

% %Down-sample params
% dwn_smpl_F = 250;
% dwn_smpl_T = 1/dwn_smpl_F;
% dwn_smpl = Fs/dwn_smpl_F;

%% PAC signal params

SNR = 3;
n_smpls = 100;
rng(200,"twister"); rng_seeds = randi([1,1e5],n_smpls,1);

% PAC model params
LF_freq_1 = [6,7]; HF_freq_1 = [55,60];
disp(['LF = ', num2str(LF_freq_1), ', HF1 = ', num2str(HF_freq_1)]);
m = 0.5; A = 200; C = 1*1e-6; am_lag = 16;

% PAC signal trimming params
tm_windw = 10; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
len_trim_ind = length(trim_ind);

%%% PAC
s_final_mat = zeros(len_trim_ind, n_smpls);

parfor i = 1:n_smpls
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

    %s_final = fft_bndpss_flt( s_final, Fs, 0, dwn_smpl_F/2 );
    %s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
    %s_final = lowpass_fir(s_final, (dwn_smpl_F/2)-2, Fs);
    %s_final = s_final(1:dwn_smpl:N);
    %tspan = tspan(1:dwn_smpl:N);
    %Fs = dwn_smpl_F; Ts = 1/Fs;

    %pink = pink(1:dwn_smpl:N);
    %=======================

    %=======================
    %Trim the PAC signal to get a small segment

    s_final_trim = s_final(trim_ind)';
    %N = length(s_final_trim);
    %----------------------
    fftn = length(s_final_trim);
    %tm_frq_plt(s_final_trim, Fs, fftn);
    %=======================

    %=======================
    % Add pink noise

    pink         = pink(trim_ind)';
    s_final_trim = ( SNR*( std(pink)/std(s_final_trim) ) ) .* s_final_trim;
    s_final_trim = pink + s_final_trim;
    %=======================

    s_final_mat(:,i) = s_final_trim;
    %% 

end

%% ------------------------ Filtering-based PAC Comodulogram ------------------------------

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
%%

%{
figure; tiledlayout(2,2);
nexttile;
imagesc(fL_vals, fH_vals, mean(Comods_mat_OzktMI,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
nexttile;
imagesc(fL_vals, fH_vals, mean(Comods_mat_CanltyMVL,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
nexttile;
imagesc(fL_vals, fH_vals, mean(Comods_mat_TortMI,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
nexttile;
imagesc(flow_GLM, fhigh_GLM, mean(Comods_mat_GLM_robust,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
%}

figure; tiledlayout(2,2);
nexttile;
imagesc(fL_vals, fH_vals, median(Comods_mat_OzktMI,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);

nexttile;
imagesc(fL_vals, fH_vals, median(Comods_mat_CanltyMVL,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);

nexttile;
imagesc(fL_vals, fH_vals, median(Comods_mat_TortMI,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);

nexttile;
imagesc(flow_GLM, fhigh_GLM, median(Comods_mat_GLM_robust,3)'); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rect_pos_1 = [LF_freq_1(1)-0.25, HF_freq_1(1)-0.5, (abs(diff(LF_freq_1))*0.5)+1, (abs(diff(HF_freq_1))*1)+1]; rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);


%%
%%

%% Local functions - PAC general

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



