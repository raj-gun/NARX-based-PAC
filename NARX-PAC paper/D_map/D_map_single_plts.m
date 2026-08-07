%D_MAP_SINGLE_PLTS Assemble the three discriminator maps shown in Figure 17.
%   Saved results for the 7/63 Hz, 6-7/55-60 Hz, and two-band experiments
%   are loaded into separate structures. The script extracts each sign map
%   D and plots them on a common horizontal layout.
%
clear;close all;clc;
%% Load three MAT-file workspaces into separate structure variables
NARX_PAC_7_63       = load('\<path-to>\NARX-PAC paper\7_PAC_63\Fig 13\7-63_pinknoise.mat');
NARX_PAC_rng_single = load('\<path-to>\NARX-PAC paper\6-7_PAC_55-65\6-7_55-65_pinknoise_SNR_3.mat');
NARX_PAC_rng_double = load('\<path-to>\NARX-PAC paper\9-10_PAC_35-40_n_70-80\9-10_35-40_70-80_pinknoise_SNR_3.mat');


%% Extract x, y, and z arrays
% Change these field names to match the variables stored in each MAT-file.
x1 = NARX_PAC_7_63.fL_vals;
y1 = NARX_PAC_7_63.fH_vals;
z1 = NARX_PAC_7_63.comod_D;

x2 = NARX_PAC_rng_single.fL_vals;
y2 = NARX_PAC_rng_single.fH_vals;
z2 = NARX_PAC_rng_single.comod_D;

x3 = NARX_PAC_rng_double.fL_vals;
y3 = NARX_PAC_rng_double.fH_vals;
z3 = NARX_PAC_rng_double.comod_D;

%% Create three horizontally arranged imagesc plots
figure;

tiledlayout(1, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

font_size = 25;

nexttile;
imagesc(x1, y1, z1);
set(gca, 'YDir', 'normal');
axis tight;
colorbar;
xlabel('Low Frequency (Hz)');
ylabel('High Frequency (Hz)');
set(gca, 'FontSize', font_size);

nexttile;
imagesc(x2, y2, z2);
set(gca, 'YDir', 'normal');
axis tight;
colorbar;
xlabel('Low Frequency (Hz)');
set(gca, 'FontSize', font_size);

nexttile;
imagesc(x3, y3, z3);
set(gca, 'YDir', 'normal');
axis tight;
colorbar;
xlabel('Low Frequency (Hz)');
set(gca, 'FontSize', font_size);

