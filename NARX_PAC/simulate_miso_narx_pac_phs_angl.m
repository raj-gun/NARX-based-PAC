%SIMULATE_MISO_NARX_PAC_PHS_ANGL Build preferred-phase maps from saved models.
%   This post-identification script loads NARX-PAC results, decomposes every
%   retained model into Sigma_u1 and Sigma_u2+Sigma_u1u2, computes the
%   phase-amplitude histogram for each frequency pair, and groups the results
%   by low-frequency centre. The resulting maps support the preferred-phase
%   characterisation described in Section IV and Figures 14 and 19.
%   Required loaded variables include Comods, fL_vals, fH_vals, fL_grd,
%   fH_grd, s_final, s_final_trim, Fs, Ts, N, frq_bndw_LF,
%   frq_bndw_HF, filt_typ, All_freq_comb, and narx_pac_modls_2.
%   This file is a script, not a function. Update the absolute paths before
%   use. It expects previously identified models and does not run the complete
%   NARX-PAC grid search.

%% Reset the workspace and configure paths
clear all;clc;close all;


addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\NonSysID-i\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\Utils');
%% Load and inspect the saved PAC result
file_dir = 'C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\Results\MISO_NARX\Tortlab\';
file_name = '500Hz_200s-0wndw_1-0.5-sbp_sbp_LFcos_2';
save_file_name = 'phs_PltDat_2.mat';
save_file_dir = file_dir;%'C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\Results\MISO_NARX\Tortlab\Figs\';
save_plt = false;
load([file_dir,file_name,'.mat']);
Comod = Comods{1};
figure; imagesc(fL_vals, fH_vals, Comod); colorbar; axis xy; set(gca, 'FontSize', 18); %axis([3,20, 35,95,-inf,inf]);
figure; surf(fL_grd, fH_grd, Comod); shading interp; %axis([2,20, 25.5,90,-inf,inf]);
tm_frq_plt(s_final, Fs, length(s_final));
%% Define numerical helpers and phase-analysis settings
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);
freq_pos = @(freq,Fs,fftn) floor(freq*fftn/Fs) + 1;

%% Phase/delay analysis
%% Simulate canonical components and compute phase-amplitude histograms
query_freqs = unique(All_freq_comb(:,1));
n_query_freqs = length(query_freqs);
fft_res = Fs/N;
K_LF = round(frq_bndw_LF/fft_res); if K_LF==0; K_LF=1; end %floor( ( ( (frq_bndw_LF)/Fs )*N ) );
K_HF = round(frq_bndw_HF/fft_res); if K_HF==0; K_HF=1; end %floor( ( ( (frq_bndw_HF)/Fs )*N ) );
nbins = 100;
phi = [0,0];
plot_data = cell(n_query_freqs,1);

cnt=1;
for query_freq = query_freqs'
    ind_phs_ana = find(All_freq_comb(:,1)==query_freq);
    n_itms = length(ind_phs_ana);
    amp_mean_mat = zeros(nbins,n_itms);
    bin_centers_mat = zeros(nbins,n_itms);

    for i=1:n_itms
        j = ind_phs_ana(i);
        pos_freq_comp = All_freq_comb(j,[1,2]);
        model = narx_pac_modls_2{j,1};
        %------------------------- Phase/delay analysis ------------------------------------------
        sig_filt_LF = nrrw_bnd_fft_filt(s_final_trim', Fs, pos_freq_comp(1), K_LF, filt_typ{1});
        sig_filt_HF = nrrw_bnd_fft_filt(s_final_trim', Fs, pos_freq_comp(2), K_HF, filt_typ{2});

        std_cos = sqrt(0.5);
        sig_filt_std = std([sig_filt_LF', sig_filt_HF'] ,1); scl_fctr = sig_filt_std./std_cos;

        [LF_sig, AM_sig] = model_simulation_OSD(model,pos_freq_comp,Ts,phi,scl_fctr);
        LF_phs_mod = angle(hilbert(LF_sig));
        HF_env_mod = abs(hilbert(AM_sig));
        [bin_centers, amp_means] = pac_hist(nbins , LF_phs_mod , HF_env_mod , false);
        amp_mean_mat(:,i) = amp_means; bin_centers_mat(:,i) = bin_centers;

        if size(LF_phs_mod,2) ~= 1; LF_phs_mod = LF_phs_mod.'; end % Transpose to a column vector
        if size(HF_env_mod,2) ~= 1; HF_env_mod = HF_env_mod.'; end % Transpose to a column vector

    end

    HF_freqs = All_freq_comb(ind_phs_ana,2);
    [phs_grd, HF_freqs_grd] = meshgrid(bin_centers_mat(:,1), HF_freqs);


    plot_data{cnt,1} = {bin_centers_mat, HF_freqs, amp_mean_mat, phs_grd , HF_freqs_grd};

    cnt=cnt+1;
end

%% Package and optionally save the preferred-phase data
plot_data = {plot_data , query_freqs};
if save_plt
    save([save_file_dir, save_file_name], 'plot_data');
end

%% Original signal PAC phase analysis


%% Visualise preferred-phase structure for each low-frequency query
cnt = 1;
for query_freq = query_freqs'

    figure; sgtitle(num2str(query_freq));
    subplot(1,2,1);
    imagesc(bin_centers_mat(:,1), HF_freqs, amp_mean_mat'); axis xy;
    set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    set(gca, 'FontSize', 20);

    subplot(1,n_query_freqs,cnt); title(num2str(query_freq));
    imagesc(bin_centers_mat(:,1), HF_freqs, amp_mean_mat'); axis xy;
    set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    set(gca, 'FontSize', 20);

    if n_itms > 1
        subplot(1,2,2);
        surf(phs_grd , HF_freqs_grd , amp_mean_mat'); shading interp;
        set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        axis tight
        set(gca, 'FontSize', 20);
    end

    cnt = cnt + 1;
end


%% Local functions
function [bin_centers, amp_means] = pac_hist(nbins,LF_phs,HF_env, plt)
%PAC_HIST Compute mean high-frequency amplitude within slow-phase bins.
%   Inputs are the number of bins, slow instantaneous phase in radians,
%   high-frequency amplitude envelope, and an optional plotting flag.
%   Outputs are phase-bin centres and the corresponding mean amplitudes.
edges = linspace(-pi, pi, nbins+1); % Phase bins
amp_means = zeros(1, nbins);        % To store mean amplitudes
% Bin the phase data and compute mean amplitude in each bin
for k = 1:nbins
    indices = LF_phs >= edges(k) & LF_phs < edges(k+1);
    amp_means(k) = mean(HF_env(indices));
end
% For visualization: bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
if plt
    figure;
    bar(rad2deg(bin_centers), amp_means, 'FaceColor', [0.2 0.6 0.8], 'EdgeAlpha', 0);
    set(gca,"TickLabelInterpreter",'latex'); set(gca, 'XTick', [-180 -90 0 90 180]); set(gca, 'XTickLabel', {'$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'});
    set(gca,'YTick',[]);
    set(gca, 'FontSize', 20);
end
end
