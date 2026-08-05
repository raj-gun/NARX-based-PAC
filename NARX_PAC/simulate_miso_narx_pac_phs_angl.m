clear all;clc;close all;
% addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\NonSysID\NonSysID-main\NonSysID\');
% addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\NonSysID-i\NonSysID-i\');
% addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NonSysID_cpy23_i\');
% addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NonSysID_cpy23_i\Utils');


addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\NonSysID-i\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\NARX_PAC\Utils');
%%
% [file_name,file_dir] = uigetfile;
file_dir = 'C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\Results\MISO_NARX\Tortlab\';
% file_dir = 'C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\Results\MISO_NARX\Tortlab\Fixed bandwidth\Delay_0o5\Lin_PRESS_3_NL_PRESS_6\Aft Chngng Condns\';
file_name = '500Hz_200s-0wndw_1-0.5-sbp_sbp_LFcos_2';
%-----------------
save_file_name = 'phs_PltDat_2.mat';
save_file_dir = file_dir;%'C:\Users\rajin\OneDrive - Coventry University\PhD project\Matlab files\CFC\Results\MISO_NARX\Tortlab\Figs\';
save_plt = false;
%-----------------
load([file_dir,file_name,'.mat']);
Comod = Comods{1};
figure; imagesc(fL_vals, fH_vals, Comod); colorbar; axis xy; set(gca, 'FontSize', 18); %axis([3,20, 35,95,-inf,inf]);
figure; surf(fL_grd, fH_grd, Comod); shading interp; %axis([2,20, 25.5,90,-inf,inf]);
tm_frq_plt(s_final, Fs, length(s_final));
%%
approx = @(value,acc) round(value/acc)*acc;
round_up = @(value,acc) floor(value) + ceil( (value-floor(value))/acc) * acc;
rand_rng = @(a,b) a + (b-a)*rand;
rand_rng_arry = @(a,b,c) a + (b-a)*rand(c,1);
freq_pos = @(freq,Fs,fftn) floor(freq*fftn/Fs) + 1;

%% Phase/delay analysis
query_freqs = unique(All_freq_comb(:,1));
n_query_freqs = length(query_freqs);
fft_res = Fs/N;
K_LF = round(frq_bndw_LF/fft_res); if K_LF==0; K_LF=1; end %floor( ( ( (frq_bndw_LF)/Fs )*N ) );
K_HF = round(frq_bndw_HF/fft_res); if K_HF==0; K_HF=1; end %floor( ( ( (frq_bndw_HF)/Fs )*N ) );
nbins = 100;
phi = [0,0];
plot_data = cell(n_query_freqs,1);
% phs_angles = angle(fft(s_final_trim,N));

% figure;
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
        %----------------------------------------
        % phi = phs_angles(floor(pos_freq_comp./fft_res) + 1);
        %------------------------- Phase/delay analysis ------------------------------------------
        sig_filt_LF = nrrw_bnd_fft_filt(s_final_trim', Fs, pos_freq_comp(1), K_LF, filt_typ{1});
        sig_filt_HF = nrrw_bnd_fft_filt(s_final_trim', Fs, pos_freq_comp(2), K_HF, filt_typ{2});

        std_cos = sqrt(0.5);
        sig_filt_std = std([sig_filt_LF', sig_filt_HF'] ,1); scl_fctr = sig_filt_std./std_cos;

        %scl_fctr(1) = 1;
        % scl_fctr = [1,1];
        [LF_sig, AM_sig] = model_simulation_OSD(model,pos_freq_comp,Ts,phi,scl_fctr);
        %----------------------------------------
        LF_phs_mod = angle(hilbert(LF_sig));
        HF_env_mod = abs(hilbert(AM_sig));
        [bin_centers, amp_means] = pac_hist(nbins , LF_phs_mod , HF_env_mod , false);
        amp_mean_mat(:,i) = amp_means; bin_centers_mat(:,i) = bin_centers;

        if size(LF_phs_mod,2) ~= 1; LF_phs_mod = LF_phs_mod.'; end % Transpose to a column vector
        if size(HF_env_mod,2) ~= 1; HF_env_mod = HF_env_mod.'; end % Transpose to a column vector

        %LF_phs_mod_mat(:,i) = LF_phs_mod;
        %HF_env_mod_mat(:,i) = HF_env_mod;
        % ----------------------------------------
        %-----------------------------------------
    end

    HF_freqs = All_freq_comb(ind_phs_ana,2);
    [phs_grd, HF_freqs_grd] = meshgrid(bin_centers_mat(:,1), HF_freqs);

    % figure; sgtitle(num2str(query_freq));
    % subplot(1,2,1);
    % imagesc(bin_centers_mat(:,1), HF_freqs, amp_mean_mat'); axis xy;
    % set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    % set(gca, 'FontSize', 20);
    % 
    % subplot(1,n_query_freqs,cnt); title(num2str(query_freq));
    % imagesc(bin_centers_mat(:,1), HF_freqs, amp_mean_mat'); axis xy;
    % set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    % set(gca, 'FontSize', 20);
    % 
    % if n_itms > 1
    %     % subplot(1,2,2);
    %     figure; sgtitle(num2str(query_freq));
    %     surf(bin_centers_mat' , HF_freqs_grd , amp_mean_mat','EdgeColor','none'); %shading interp;
    %     set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    %     axis tight
    %     set(gca, 'FontSize', 20); view(2);
    % end

    plot_data{cnt,1} = {bin_centers_mat, HF_freqs, amp_mean_mat, phs_grd , HF_freqs_grd};

    cnt=cnt+1;
end

plot_data = {plot_data , query_freqs};
if save_plt
    save([save_file_dir, save_file_name], 'plot_data');
end

%% Original signal PAC phase analysis
%{
LF_phs_mod = angle(hilbert(s_LF));
HF_env_mod = abs(hilbert(s_HF1_shft));
[bin_centers, amp_means] = pac_hist(nbins , LF_phs_mod , HF_env_mod , true);
%% Plot original data without noise and NARX-based PAC
freq_ind = 2;
figure;
bar(rad2deg(plot_data{1,1}{freq_ind,1}{1,1}) , plot_data{1,1}{freq_ind,1}{1,3}, 'FaceColor', [0.47,0.67,0.19], 'EdgeAlpha', 0);
set(gca,"TickLabelInterpreter",'latex'); set(gca, 'XTick', [-180 -90 0 90 180]); set(gca, 'XTickLabel', {'$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'});
set(gca,'YTick',[]);
set(gca, 'FontSize', 20);
%}

% plot( plot_data{1,1}{freq_ind,1}{1,1} , plot_data{1,1}{freq_ind,1}{1,3}); hold on
% plot(bin_centers,amp_means, 'LineWidth', 2); hold off;
% legend( [ num2cell( num2str( plot_data{1,1}{freq_ind,1}{1,2} ) ,2) ;  'Original without noise'] );
% title( num2str( plot_data{1,2}(freq_ind) ) )

%%
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
        % figure; sgtitle(num2str(query_freq));
        surf(phs_grd , HF_freqs_grd , amp_mean_mat'); shading interp;
        set(gca, 'XTick', [-pi, -0.5*pi, 0, 0.5*pi, pi]); set(gca, 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        axis tight
        set(gca, 'FontSize', 20);
    end

    cnt = cnt + 1;
end

%%


%%
%% Local functions
%%
function [bin_centers, amp_means] = pac_hist(nbins,LF_phs,HF_env, plt)
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
    % xlabel('Phase of low frequency (rads)');
    % ylabel('Amplitude of high frequency');
    set(gca,"TickLabelInterpreter",'latex'); set(gca, 'XTick', [-180 -90 0 90 180]); set(gca, 'XTickLabel', {'$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'});
    set(gca,'YTick',[]);
    set(gca, 'FontSize', 20);
end
end