clear all;clc;close all;

addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\NonSysID-i\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\Utils\');
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

%% Non-sinusoidal

LF = 7;
u = sin( (2*pi*LF).*tspan )./LF;
s_final = 10 ./ ( 1 + exp( -12.*( 0.5.*(1+5.*u) -0.7 ) )  );
am_lag = 0;

% s_final = van_d_pol_LF(tspan,116.5)';
% am_lag = 0;

n_smpls = 100; rng(100,"twister"); rng_seeds = randi([1,1e5],n_smpls,1); 
rng( rng_seeds(60) ); [~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, [9.5,10.5], [55,65]);
% s_final = 0.041.*s_final;
% sn_ratio = snr(s_final,pink); disp(['SNR = ', num2str(sn_ratio)]);
% s_final = pink + s_final;
% tm_frq_plt( s_final, Fs, N);

%% Down sample

dwn_smpl_F = 250;
%s_final = fft_bndpss_flt( s_final, Fs, 0, dwn_smpl_F/2 );
s_final = lowpss_fft_filt(s_final, Fs, (dwn_smpl_F/2)-2, 2);
% s_final = lowpass_fir(s_final, dwn_smpl_F/2, Fs);
dwn_smpl = Fs/dwn_smpl_F;
s_final = s_final(1:dwn_smpl:N);
tspan = tspan(1:dwn_smpl:N);
Fs = dwn_smpl_F; Ts = 1/Fs;

pink = pink(1:dwn_smpl:N);

% N = length(tspan);
% fftn = 1000;%Fs/N;
% w = 0:Fs/fftn:Fs-(Fs/fftn);

%% Test single  sample of noise

%Trim the PAC signal to get a small segment
tm_windw = 20; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
s_final_trim = s_final(trim_ind)';%(600:1100)';

%{1
rng(100,"twister");
SNR = 3;
pink = pink(trim_ind)';
s_final_trim = ( SNR*( std(pink)/std(s_final_trim) ) ) .* s_final_trim;
sn_ratio = snr(s_final_trim,pink); disp(['SNR = ', num2str(sn_ratio), 'dB or ',  num2str(db2mag(sn_ratio))]);
s_final_trim = pink + s_final_trim;
%}

N = length(s_final_trim);
figure;plot(tspan(trim_ind),s_final_trim);
%----------------------
fftn = length(s_final_trim); 
w = 0:Fs/fftn:Fs-(Fs/fftn);
% s_final_trim_fft = Ts.*fft(s_final_trim, fftn);figure;subplot(2,1,1);plot(w, abs(s_final_trim_fft));subplot(2,1,2);plot(w, angle(s_final_trim_fft).*(180/pi));
tm_frq_plt(s_final_trim, Fs, fftn);

fL_vals = [1:1:20]; fH_vals = [15:1:90];
% fL_vals = [1:0.5:15]; fH_vals = [18:1:50];
% fL_vals = [5:0.5:20]; fH_vals = [30:1:100];


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

% [HF_MI_score_1, HF_MI_score_final] = SpuCup_intrmd(fL_vals, fH_vals, diff_comod);

[Comod_harmonic_rmv, IF_harmonic_test_dat] = IF_harmonic_test(All_freq_comb_1, phs_data_mat, fL_vals, fH_vals, Comods{1}, Ts);
figure; imagesc(fL_vals, fH_vals, Comod_harmonic_rmv); colorbar; axis xy; set(gca, 'FontSize', 18);
sgtitle('Commod after removing harmonics');

%%
figure; imagesc(fL_vals, fH_vals, Comods{1}); colorbar; axis xy; set(gca, 'FontSize', 18);

% figure; imagesc(fL_vals, fH_vals, diff_comod); colorbar; axis xy; set(gca, 'FontSize', 18);

% figure; surf(fL_grd, fH_grd, Comod, 'EdgeColor','none'); axis tight; view(2); axis tight; shading interp;

% disp(All_freq_comb( All_freq_comb(:, 1)==fL & All_freq_comb(:, 2)==fH , :))
%%
%%

%% =====================================================
%% Local functions
%% =====================================================

%% Non-sinusoidal signal
function [s_LF] = van_d_pol_LF(tspan,w0)
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