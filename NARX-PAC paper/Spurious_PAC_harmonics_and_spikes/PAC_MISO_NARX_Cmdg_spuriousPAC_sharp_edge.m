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
%% Spurious PAC
%% =====================================================

%%
tspan = 0:Ts:25-Ts;%(N*Ts-Ts);
N = length(tspan);
fftn = 4000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);

%% Sharp edge
LF = 9;
s_final = sharp_edge(LF, 0.2, Ts, N);

n_smpls = 100; rng(100,"twister"); rng_seeds = randi([1,1e5],n_smpls,1); 
rng( rng_seeds(60) ); [~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, [9.5,10.5], [55,65]);

%% Down sample

dwn_smpl_F = 250;
s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
dwn_smpl = Fs/dwn_smpl_F;
s_final = s_final(1:dwn_smpl:N);
tspan = tspan(1:dwn_smpl:N);
Fs = dwn_smpl_F; Ts = 1/Fs;

pink = pink(1:dwn_smpl:N);


%% Test single  sample of noise

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


%% ------------------------ NARX based MISO PAC Comodulogram ------------------------------
filt_typ = {'sbp','sbp'}; % bw , sbp , guss
frq_bndw_LF = 0.5;
frq_bndw_HF = 0.5;
disp(['frq_bndw_LF = ', num2str(frq_bndw_LF), ', frq_bndw_HF = ', num2str(frq_bndw_HF)]);

tic
[Comods , diff_comod, phs_data_mat, fL_grd, fH_grd, All_freq_comb_1, All_freq_comb, All_freq_comb_ARX_1, All_freq_comb_ARX_2, narx_pac_modls_1 , narx_pac_modls_2]...
    = pac_miso_Cmdg_mod_21(s_final_trim, fL_vals, fH_vals, Fs, 3, filt_typ, frq_bndw_LF, frq_bndw_HF);
toc
%% Post processing of results


[Comod_harmonic_rmv, IF_harmonic_test_dat] = IF_harmonic_test(All_freq_comb_1, phs_data_mat, fL_vals, fH_vals, Comods{1}, Ts);
figure; imagesc(fL_vals, fH_vals, Comod_harmonic_rmv); colorbar; axis xy; set(gca, 'FontSize', 18);
sgtitle('Commod after removing harmonics');

%%
figure; imagesc(fL_vals, fH_vals, Comods{1}); colorbar; axis xy; set(gca, 'FontSize', 18);



%%
%%

%% =====================================================
%% Local functions
%% =====================================================

%% Sharp edges
% Code is adapted from Kramer et al. (2008), Jrn. Nrsc. Methds. 
% and Ozkurt et al., (2011) Jrn. Nrsc. Methds.
function [s] = sharp_edge(f, edge_pos, Ts, N)

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
