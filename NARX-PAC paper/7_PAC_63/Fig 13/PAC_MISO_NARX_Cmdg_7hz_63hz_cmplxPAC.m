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
%% PAC LF-sine HF-sine simple model

tspan = 0:Ts:25-Ts;%(N*Ts-Ts);
N = length(tspan);
fftn = 1000;%Fs/N;
w = 0:Fs/fftn:Fs-(Fs/fftn);

fL=7; fH=63;
f_phi = [39.4212  290.5146];
am_lag = 60;
%=======================
s_LF = cos((2*pi*fL).*tspan + (f_phi(1)*pi/180));
s_HF = cos((2*pi*fH).*tspan + (f_phi(2)*pi/180));

d_F = 0.2;
LF_freq = [-d_F,d_F]+fL; HF_freq = ([-d_F,d_F])+fH;
rng_i = 1000;
rng(rng_i+1000); [~, ~, pink, ~] = pink_noise_LF_HF(N, Fs, LF_freq, HF_freq);
%=======================
% Non-sine PAC
%--------
m = 0.5;
[s_final, s_HF1_shft] = pac_general_1(s_LF, s_HF, 200, 1*1e-6, m, am_lag);
%--------
%=======================
disp([fL,fH]); disp([f_phi,am_lag]);

tm_frq_plt(s_final, Fs, N);
%% Visualise PAC signal
if am_lag~=0
    s_final = s_final(am_lag:end);
    s_LF = s_LF(am_lag:end);
    s_HF1_shft = s_HF1_shft(am_lag:end);
    tspan = tspan(am_lag:end);
    N = length(tspan);
    fftn = 1000;%Fs/N;
    w = 0:Fs/fftn:Fs-(Fs/fftn);
end
s_final_org = s_final;
s_final_org_fft = Ts.*fft(s_final_org, fftn);
figure;subplot(2,1,1);plot(w, abs(s_final_org_fft));subplot(2,1,2);plot(w, angle(s_final_org_fft).*(180/pi));
%% Down sample

dwn_smpl_F = 250;
s_final = lowpass_fir(s_final, (dwn_smpl_F/2)-2, Fs);
dwn_smpl = Fs/dwn_smpl_F;
s_final = s_final(1:dwn_smpl:N);
tspan = tspan(1:dwn_smpl:N);
s_LF = s_LF(1:dwn_smpl:N);
s_HF1_shft = s_HF1_shft(1:dwn_smpl:N);
Fs = dwn_smpl_F; Ts = 1/Fs;

pink = pink(1:dwn_smpl:N);


%% Test single  sample of noise

%Trim the PAC signal to get a small segment
tm_windw = 10; tm_itr = 0;
trim_ind = (tm_itr*tm_windw/Ts)+1:((tm_itr+1)*tm_windw)/Ts;%(1 + (30/fL)/Ts);%
s_final_trim = s_final(trim_ind)';%(600:1100)';

%% Add pink noise
pink = pink(trim_ind)';
s_final_trim = ( 3*( std(pink)/std(s_final_trim) ) ) .* s_final_trim;
sn_ratio = snr(s_final_trim,pink); disp(['SNR = ', num2str(sn_ratio), 'dB or ',  num2str(db2mag(sn_ratio))]);
tm_frq_plt(s_final_trim, Fs, fftn);
s_final_trim = pink + s_final_trim;
tm_frq_plt(s_final_trim, Fs, fftn);
%%
N = length(s_final_trim);
figure;plot(tspan(trim_ind),s_final_trim);
%----------------------
fftn = length(s_final_trim); 
w = 0:Fs/fftn:Fs-(Fs/fftn);
tm_frq_plt(s_final_trim, Fs, fftn);

fL_vals = [4:1:10]; fH_vals = [30:1:100];


%% ------------------------ NARX based MISO PAC Comodulogram ------------------------------
filt_typ = {'sbp','sbp'}; % bw , sbp , guss
frq_bndw_LF = 0.25; 
frq_bndw_HF = 0.5; %35-35
disp(['frq_bndw_LF = ', num2str(frq_bndw_LF), ', frq_bndw_HF = ', num2str(frq_bndw_HF)]);

tic
[Comods , diff_comod, phs_data_mat, fL_grd, fH_grd, All_freq_comb_1, All_freq_comb, All_freq_comb_ARX_1, All_freq_comb_ARX_2, narx_pac_modls_1 , narx_pac_modls_2]...
    = pac_miso_Cmdg_mod_21(s_final_trim, fL_vals, fH_vals, Fs, 3, filt_typ, frq_bndw_LF, frq_bndw_HF);
toc

%%
fL_diff = mean(abs(diff(fL_vals))); fH_diff = mean(abs(diff(fH_vals))); 
rect_pos_box = @(LF_freq, HF_freq, fL_diff, fH_diff) [LF_freq(1)-fL_diff*0.5, HF_freq(1)-fH_diff*0.5, (abs(diff(LF_freq))*fL_diff)+1, (abs(diff(HF_freq))*1)+1]; 

figure; imagesc(fL_vals, fH_vals, Comods{1}); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;

comod_D = diff_comod;
comod_D( comod_D > 0 ) = 1; comod_D( comod_D < 0 ) = -1; 
figure; imagesc(fL_vals, fH_vals, comod_D); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;

%%  Post processing of results


[Comod_harmonic_rmv, IF_harmonic_test_dat] = IF_harmonic_test (All_freq_comb_1, phs_data_mat, fL_vals, fH_vals, Comods{1}, Ts);
figure; imagesc(fL_vals, fH_vals, Comod_harmonic_rmv); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;

[HF_MI_score_1_2, HF_MI_score_final_2, Comod_intrmd] = SpuCup_intrmd_2(fL_vals, fH_vals, Comods{1}, diff_comod, All_freq_comb_1, phs_data_mat, Fs, 1);
figure; imagesc(fL_vals, fH_vals, Comod_intrmd); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;

[Comod_harmonic] = AM_var_commod(All_freq_comb_1, phs_data_mat, fL_vals, fH_vals);
figure; imagesc(fL_vals, fH_vals, Comod_harmonic); colorbar; axis xy; set(gca, 'FontSize', 18);


%% Functions Used
%% Creat LF and HF signals

function [s_LF, s_HF] = rand_varying_LF_HF(LF_freq, HF_freq, Fs, N)
s_LF = randn(1,N);
s_LF = fft_bndpss_flt( s_LF, Fs, LF_freq(1), LF_freq(2) ); s_LF = s_LF./max(abs(s_LF));

s_HF = randn(1,N);
s_HF = fft_bndpss_flt( s_HF, Fs, HF_freq(1), HF_freq(2) ); s_HF = ( s_HF./max(abs(s_HF)) ).*0.25;
end

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

%% Other functions

%FFT base bandpass filter f1<f2
function data = fft_bndpss_flt(data_raw,Fs,f1,f2)

data = bandpass(data_raw , [f1,f2] , Fs);
end

function [F_phs] = xcorr_phs_estm(F, s_final_trim, Ts)
tspan_corr = 0:Ts:(length(s_final_trim)*Ts-Ts);%(1/F)+Ts;%
phi_prb = (-pi:0.01:pi)'; phi_pprb_len = length(phi_prb);
prb_inpt = cos( ( (2*pi).*F.*tspan_corr ) + phi_prb ); % probing input

xcorr_val = zeros(phi_pprb_len,1);

for i=1:phi_pprb_len; [xcorr_val_tmp , ~] = xcorr(s_final_trim, prb_inpt(i,:) ,0); xcorr_val(i) = xcorr_val_tmp; end
[max_corr_phi_prb , max_corr_ind_phi_prb] = max(xcorr_val);
F_phs = phi_prb(max_corr_ind_phi_prb);
end

function [data_noise] = add_noise(noise_type, snr, data)
switch noise_type
    case 'white'
        %White noise
        wn = randn(size(data));
        wn = wn ./ sqrt( snr*(sum(wn.^2) / sum(data.^2)) );
        data_noise = data + wn;
    case 'pink'
        %Pink noise
        pn = pinknoise( length(data) , 1);
        pn = pn ./ sqrt( snr*(sum(pn.^2) / sum(data.^2)) );
        data_noise = data + pn;
end
end
%% PAC general

% Code adapted from Jiang et al., (2015) NeuroImage
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
%% Surrogate

function [surr] = randm_swap_surrg(data, nsurr, corr_min)
dat_len = length(data);
surr = zeros(dat_len,nsurr);
parfor i=1:nsurr
    while 1
        rng shuffle;
        surr_dat = data( randperm(dat_len) , 1 );
        corr_mat = corrcoef([surr_dat, data]);
        corr_surr_dat = abs(corr_mat(1,2));
        if corr_surr_dat <= corr_min
            surr(:,i) = surr_dat;
            break;
        end
    end
end
end
