
% 
% file_dir_PP = '\<path-to>\Results\MISO_NARX\SyntheticData\Multisample\';
% file_name_PP = 'F1_3_5s_wrk.mat';
% load([file_dir_PP,file_name_PP]);

%%
Comods = cellfun(@(x) x(1), Comods);

Comods_intrmd = cell(1,n_smpls);
Comods_harm = cell(1,n_smpls);
Comods_D = cell(1,n_smpls);
for i = 1:n_smpls
    phs_data_mat = PAC_data_mat{1,i}; 
    All_freq_comb_1 = All_freq_comb_dat{1,i};
    Comod = Comods{1,i};
    diff_comod = Comods_diff{1,i};
    if isempty(All_freq_comb_1);disp(i);  continue; end

    [~, ~, Comod_intrmd_rmv] = SpuCup_intrmd_2(fL_vals, fH_vals, Comod, diff_comod, All_freq_comb_1, phs_data_mat, Fs, 0);

    [Comod_harmonic_rmv, ~] = IF_harmonic_test(All_freq_comb_1, phs_data_mat, fL_vals, fH_vals, Comod, Ts);

    comod_D = diff_comod;
    comod_D( comod_D > 0 ) = 1; comod_D( comod_D < 0 ) = -1; 

    Comods_intrmd{1,i}  = Comod_intrmd_rmv;
    Comods_harm{1,i}    = Comod_harmonic_rmv;
    Comods_D{1,i}       = comod_D;

end

%%

Comods_mat          = cat(3, Comods{:});
Comods_intrmd_mat   = cat(3, Comods_intrmd{:});
Comods_harm_mat     = cat(3, Comods_harm{:});
Comods_diff_mat     = cat(3, Comods_diff{:});
Comods_D_mat        = cat(3, Comods_D{:});

%%

% file_dir_PP = '\<path-to>\Results\MISO_NARX\SyntheticData\Multisample\PP\';
% file_name_PP = [ file_name_PP(1:end-4) , '_PP' , '.mat' ];
% save([file_dir_PP,file_name_PP]);
%%

fL_diff = mean(abs(diff(fL_vals))); fH_diff = mean(abs(diff(fH_vals))); 
rect_pos_box = @(LF_freq, HF_freq, fL_diff, fH_diff) [LF_freq(1)-fL_diff*0.5, HF_freq(1)-fH_diff*0.5, (abs(diff(LF_freq))*1)+1, (abs(diff(HF_freq))*1)+1]; 
rect_pos_1 = rect_pos_box(LF_freq_1, HF_freq_1, fL_diff, fH_diff);

%====================================
figure;
subplot(2,2,1); imagesc(fL_vals, fH_vals, mean(Comods_mat,3)); colorbar; axis xy; set(gca, 'FontSize', 18);
title('Commod');
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);

subplot(2,2,2); imagesc(fL_vals, fH_vals, mean(Comods_intrmd_mat,3)); colorbar; axis xy; set(gca, 'FontSize', 18);
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
title('Commod SC-i rmvd');

subplot(2,2,3); imagesc(fL_vals, fH_vals, mean(Comods_harm_mat,3)); colorbar; axis xy; set(gca, 'FontSize', 18);
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);

subplot(2,2,4); imagesc(fL_vals, fH_vals, mean(Comods_D_mat,3)); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
%====================================
figure;
subplot(2,2,1); imagesc(fL_vals, fH_vals, median(Comods_mat,3)); colorbar; axis xy; set(gca, 'FontSize', 18);
title('Commod');
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);

subplot(2,2,2); imagesc(fL_vals, fH_vals, median(Comods_intrmd_mat,3)); colorbar; axis xy; set(gca, 'FontSize', 18);
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);
title('Commod SC-i rmvd');

subplot(2,2,3); imagesc(fL_vals, fH_vals, median(Comods_harm_mat,3)); colorbar; axis xy; set(gca, 'FontSize', 18);
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);

subplot(2,2,4); imagesc(fL_vals, fH_vals, median(Comods_D_mat,3)); colorbar; axis xy; set(gca, 'FontSize', 18); hold on;
rectangle('Position',rect_pos_1, 'EdgeColor','r', 'LineWidth', 1);

% %====================================
% 
% 
% 
% %====================================
% 
% 
% 
% %====================================
