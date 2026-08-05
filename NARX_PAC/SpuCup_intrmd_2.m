function [HF_MI_score_1, HF_MI_score_final, Comod_intrmd_rmv] = SpuCup_intrmd_2(fL_vals, fH_vals, Comod, diff_comod, All_freq_comb_1, phs_data_mat, Fs, sine_freq)

%% MI and HF data

no_fH_vals = length(fH_vals);
N = size( phs_data_mat{1, 1} , 1);
% Commodulagram the grid
[fL_grd, fH_grd] = meshgrid(fL_vals, fH_vals);
pos_freq_comp_vec = [fL_grd(:) , fH_grd(:)];

HF_MI_info = [ pos_freq_comp_vec , ( pos_freq_comp_vec(:,2) + [-pos_freq_comp_vec(:,1) , pos_freq_comp_vec(:,1) ] ) ];

% [IAmp_ana_dat, IAmp_ana_dat_all] = IAmp_info(N, Fs, All_freq_comb_1, phs_data_mat, fL_vals, fH_vals);
% IAmp_ana_dat_all(:,1:4) = HF_MI_info(:,1:4);

comod_D = diff_comod;
comod_D( comod_D > 0 ) = 1; comod_D( comod_D < 0 ) = -1; 

HF_MI_info = [ HF_MI_info , comod_D(:) ]; 
% HF_MI_info = [ HF_MI_info , diff_comod(:) ];


HF_MI_score_final = HF_MI_info; HF_MI_score_final(:,end) = 0;
%% MI and HF info per low frequency
HF_MI_score_1 = zeros( size(HF_MI_score_final,1) , 3 );
HF_MI_info_2 = cell( length(fL_vals) , 1 );
commod_intrmd_HF_clstrs = cell( length(fL_vals) , 1 );

loop_cnt = 1;
for fL_freq = fL_vals
    fL_ind_lgc = HF_MI_info(:,1) == fL_freq;

    break_loop = 0;
    %% Segmentation and scoring -- MIHF

    if sine_freq == 0

        [ipt_freq, HF_MI_info_fL_1, HF_MI_info_fL_final, break_loop] = segmt_mean(HF_MI_info, fL_ind_lgc, no_fH_vals);

        %[ipt_freq, ~, ~, ~] = segmt_mean(HF_MI_info, fL_ind_lgc, no_fH_vals);
        %[HF_MI_info_fL_1, HF_MI_info_fL_final, break_loop] = segmt_mean_2(HF_MI_info, fL_ind_lgc, no_fH_vals);

        HF_MI_info_2{loop_cnt,1} = HF_MI_info_fL_1;
        commod_intrmd_HF_clstrs{loop_cnt,1} = ipt_freq;

        %     IAmp_info_fL_Idiff = IAmp_ana_dat_all(fL_ind_lgc, [1,2,6]); IAmp_info_fL_Idiff = IAmp_info_fL_Idiff( HF_MI_info_fL_nz , : );
        %     IAmp_info_fL_Ism   = IAmp_ana_dat_all(fL_ind_lgc, [1,2,7]); IAmp_info_fL_Ism   = IAmp_info_fL_Ism( HF_MI_info_fL_nz , : );
        %
        %     [IAmp_info_fL_final_Idiff, break_loop] = segmt_mean(ipt_ind, IAmp_info_fL_Idiff, IAmp_ana_dat_all(: , [1:4,6] ), fL_ind_lgc, no_fH_vals, 0);
        %     [IAmp_info_fL_final_Ism  , break_loop] = segmt_mean(ipt_ind, IAmp_info_fL_Ism  , IAmp_ana_dat_all(: , [1:4,7] ), fL_ind_lgc, no_fH_vals, 0);
        %
        %     IAmp_info_fL_final = [ IAmp_info_fL_final_Idiff(:,1:4) , IAmp_info_fL_final_Idiff(:,5) , IAmp_info_fL_final_Ism(:,5) ];
    else
        [~, HF_MI_info_fL_final, break_loop] = segmt_mean_2(HF_MI_info, fL_ind_lgc, no_fH_vals);
    end
    %%
    if break_loop == 1
        HF_MI_score_1(fL_ind_lgc,:) = [ HF_MI_info_fL_final(:, [1,2]) , 0.*HF_MI_info_fL_final(:,end)];
        loop_cnt = loop_cnt + 1;
        continue;
    end
    %%
    %--------------------------------------------------------------------
    % Positive evidence
    pos_evd_Idiff = max( [ zeros(no_fH_vals,1), -sign( HF_MI_info_fL_final(:,5) ) .* HF_MI_info_fL_final(:,6) ] , [] , 2 );
    pos_evd_Ism = max( [ zeros(no_fH_vals,1), -sign( HF_MI_info_fL_final(:,5) ) .* HF_MI_info_fL_final(:,7) ] , [] , 2  );
    pos_evd = pos_evd_Idiff + pos_evd_Ism;

    % Absolute Negative evidence
    neg_evd_HFs = sort( unique( [ HF_MI_info_fL_final( pos_evd == 2 , 3) ; HF_MI_info_fL_final( pos_evd == 2 , 4) ] ) , 'ascend');
    neg_evd = zeros(no_fH_vals,1);
    for i=1:length(neg_evd_HFs); neg_evd( HF_MI_info_fL_final(:,2) == neg_evd_HFs(i) , 1) = -1; end
    
    % Negative evidence--Intermodulations in the edges
    neg_evd_edge = -( pos_evd_Idiff.*(HF_MI_info_fL_final(:,7)==0) + pos_evd_Ism.*( HF_MI_info_fL_final(:,6)==0) );
    
    pos_evd_edge_Ism_HFs = HF_MI_info_fL_final( (pos_evd_Idiff.*(HF_MI_info_fL_final(:,7)==0))==1 , 3);
    pos_evd_edge_Ism = zeros(no_fH_vals,1);
    for i=1:length(pos_evd_edge_Ism_HFs); pos_evd_edge_Ism( HF_MI_info_fL_final(:,2) == pos_evd_edge_Ism_HFs(i) , 1) = 1; end
    %
    pos_evd_edge_Idiff_HFs = HF_MI_info_fL_final( (pos_evd_Ism.*(HF_MI_info_fL_final(:,6)==0))==1 , 4);
    pos_evd_edge_Idiff = zeros(no_fH_vals,1);
    for i=1:length(pos_evd_edge_Idiff_HFs); pos_evd_edge_Idiff( HF_MI_info_fL_final(:,2) == pos_evd_edge_Idiff_HFs(i) , 1) = 1; end
    %
    pos_evd_edge = pos_evd_edge_Idiff + pos_evd_edge_Ism;
    %--------------------------------------------------------------------

    score = 0.5.*pos_evd + 0.*pos_evd_edge + neg_evd + 0.*neg_evd_edge;

    %{
    figure;
    subplot(6,1,1); stem(HF_MI_info_fL_final(:,2),HF_MI_info_fL_final(:,5), 'r');
    subplot(6,1,2); stem( HF_MI_info_fL_final(:,2) , 0.5.*pos_evd );
    subplot(6,1,3); stem( HF_MI_info_fL_final(:,2) , 0.*pos_evd_edge );
    subplot(6,1,4); stem( HF_MI_info_fL_final(:,2) , neg_evd );
    subplot(6,1,5); stem( HF_MI_info_fL_final(:,2) , 0.*neg_evd_edge );
    subplot(6,1,6); stem( HF_MI_info_fL_final(:,2) , score );
    sgtitle(num2str(fL_freq)) ;   
    %}

    HF_MI_score_1(fL_ind_lgc,:) = [ HF_MI_info_fL_final(:, [1,2]) , score];
    %%
    %{
    HF_MI_score_1_temp = HF_MI_info_fL_final(:,1:4); HF_MI_score_1_temp = [HF_MI_score_1_temp, zeros( size(HF_MI_score_1_temp,1) ,7) ];

    lgc_nan = ( HF_MI_info_fL_final(:,3)~=0 & HF_MI_info_fL_final(:,4)~=0 );

    % Absolute HF confirm
    HF_cnfrm_Bi = ( HF_MI_info_fL_final(:,5)<0 & (HF_MI_info_fL_final(:,6)>0 & HF_MI_info_fL_final(:,7)>0) ) & lgc_nan; % MI>1  & intermodulations HF Mags. agree
    HF_cnfrm_M  = ( HF_MI_info_fL_final(:,5)>0 & (HF_MI_info_fL_final(:,6)<0 & HF_MI_info_fL_final(:,7)<0) ) & lgc_nan; % MI<1  & intermodulations HF Mags. agree
    HF_cnfrm = HF_cnfrm_Bi | HF_cnfrm_M;
    HF_MI_score_1_temp(HF_cnfrm,5) = 1;

    % Partial HF confirm
    HFp_cnfrm_Bi_Idiff = ( HF_MI_info_fL_final(:,5)<0 & HF_MI_info_fL_final(:,6)>0  ) & (lgc_nan & ~HF_cnfrm); % MI>1  & only Diff. intermodulation HF Mag. agree
    HFp_cnfrm_Bi_Ism = ( HF_MI_info_fL_final(:,5)<0 & HF_MI_info_fL_final(:,7)>0  ) & (lgc_nan & ~HF_cnfrm); % MI>1  & only Sm. intermodulation HF Mag. agree

    HFp_cnfrm_M_Idiff  = ( HF_MI_info_fL_final(:,5)>0 & HF_MI_info_fL_final(:,6)<0 ) & (lgc_nan & ~HF_cnfrm); % MI<1  & only Diff. intermodulation HF Mag. agree
    HFp_cnfrm_M_Ism  = ( HF_MI_info_fL_final(:,5)>0 & HF_MI_info_fL_final(:,7)<0 ) & (lgc_nan & ~HF_cnfrm); % MI<1  & only Sm. intermodulation HF Mag. agree
    HFp_cnfrm = (HFp_cnfrm_Bi_Idiff | HFp_cnfrm_Bi_Ism) | (HFp_cnfrm_M_Idiff | HFp_cnfrm_M_Ism);
    HF_MI_score_1_temp(HFp_cnfrm,6) = 0.5;

    % Absolute Intrmd confirm
    abs_intrmd_Idiff_vals = HF_MI_score_1_temp(HF_cnfrm, 3);
    for i=1:length(abs_intrmd_Idiff_vals); HF_MI_score_1_temp( HF_MI_score_1_temp(:,2) == abs_intrmd_Idiff_vals(i) , 7) = -1; end
    %----
    abs_intrmd_Ism_vals = HF_MI_score_1_temp(HF_cnfrm, 4);
    for i=1:length(abs_intrmd_Ism_vals); HF_MI_score_1_temp( HF_MI_score_1_temp(:,2) == abs_intrmd_Ism_vals(i) , 7) = -1; end
    %----
    %HF_MI_score_1_temp(HF_cnfrm , 7) = 0;

    %Partial intrmd confirm
    abs_intrmd_Bi_Idiff_vals = HF_MI_score_1_temp(HFp_cnfrm_Bi_Idiff, 3);
    for i=1:length(abs_intrmd_Bi_Idiff_vals); HF_MI_score_1_temp( HF_MI_score_1_temp(:,2) == abs_intrmd_Bi_Idiff_vals(i) , 8) = -0.5; end
    %----
    abs_intrmd_Bi_Ism_vals = HF_MI_score_1_temp(HFp_cnfrm_Bi_Ism, 4);
    for i=1:length(abs_intrmd_Bi_Ism_vals); HF_MI_score_1_temp( HF_MI_score_1_temp(:,2) == abs_intrmd_Bi_Ism_vals(i) , 9) = -0.5; end

    abs_intrmd_M_Idiff_vals = HF_MI_score_1_temp(HFp_cnfrm_M_Idiff, 3);
    for i=1:length(abs_intrmd_M_Idiff_vals); HF_MI_score_1_temp( HF_MI_score_1_temp(:,2) == abs_intrmd_M_Idiff_vals(i) , 10) = -0.5; end
    %----
    abs_intrmd_M_Ism_vals = HF_MI_score_1_temp(HFp_cnfrm_M_Ism, 4);
    for i=1:length(abs_intrmd_M_Ism_vals); HF_MI_score_1_temp( HF_MI_score_1_temp(:,2) == abs_intrmd_M_Ism_vals(i) , 11) = -0.5; end

    HF_MI_score_2_temp = HF_MI_info_fL_final(:,1:4); HF_MI_score_2_temp = [HF_MI_score_2_temp, zeros( size(HF_MI_score_2_temp,1) ,1) ];
    HF_MI_score_2_temp(:,end) = sum( HF_MI_score_1_temp(:,5:end) , 2);

    HF_MI_score_final(fL_ind_lgc, end) = HF_MI_score_2_temp(:,end);
    HF_MI_score_1(fL_ind_lgc, :) = HF_MI_score_1_temp;

    %{1
    figure; subplot(1,2,1); stem( HF_MI_info_fL_final(:,2) , HF_MI_info_fL_final(:,5));
    subplot(1,2,2); stem(HF_MI_score_2_temp(:,2),HF_MI_score_2_temp(:,end) , 'k' , 'LineWidth', 1.5);
    sgtitle( num2str(fL_freq) );

    figure;
    subplot(8,1,1); stem( HF_MI_info_fL_final(:,2) , HF_MI_info_fL_final(:,5));
    subplot(8,1,2); stem(HF_MI_score_1_temp(:,2),HF_MI_score_1_temp(:,5));
    subplot(8,1,3); stem(HF_MI_score_1_temp(:,2),HF_MI_score_1_temp(:,6));
    subplot(8,1,4); stem(HF_MI_score_1_temp(:,2),HF_MI_score_1_temp(:,7));
    subplot(8,1,5); stem(HF_MI_score_1_temp(:,2),HF_MI_score_1_temp(:,8));
    subplot(8,1,6); stem(HF_MI_score_1_temp(:,2),HF_MI_score_1_temp(:,9));
    subplot(8,1,7); stem(HF_MI_score_1_temp(:,2),HF_MI_score_1_temp(:,10));
    subplot(8,1,8); stem(HF_MI_score_1_temp(:,2),HF_MI_score_1_temp(:,11));
    sgtitle( num2str(fL_freq) );
    %}

    
    %%
loop_cnt = loop_cnt + 1;
end
%%

% [fL_grd, fH_grd] = meshgrid(fL_vals, fH_vals);
% pos_freq_comp_vec = [fL_grd(:) , fH_grd(:)];
% no_probes = size(pos_freq_comp_vec,1);
%
% intrmd_freq = zeros(no_probes,1);
% for i = 1:size(HF_MI_info_fL_final,1)
%     intrmd_freq(sum(pos_freq_comp_vec == HF_MI_info_fL_final(i ,[1,2]), 2) == 2, 1) =  HF_MI_info_fL_final(i ,5) ;
% end
% intrmd_freq = HF_MI_info_fL_final(:,5);
Comod_intrmd = reshape( HF_MI_score_1(:,3) , size(fL_grd)) >= 0;
Comod_intrmd_rmv = Comod .* Comod_intrmd;
% figure; imagesc(fL_vals, fH_vals, Comod_intrmd_rmv); colorbar; axis xy; set(gca, 'FontSize', 18);

% HF_MI_info_2 = cell2mat(HF_MI_info_2);
% Comod_D = reshape( HF_MI_info_2(:,3) , size(fL_grd));
% figure; imagesc(fL_vals, fH_vals, Comod_D); colorbar; axis xy; set(gca, 'FontSize', 18);

% figure;stem( HF_MI_info_fL(:,2) , HF_MI_info_fL(:,end) );
% figure;stem( HF_MI_info_fL_final(:,2) , HF_MI_info_fL_final(:,5) );
%
% figure;stem( HF_MI_info_fL(:,2) , HF_MI_info_fL(:,end) );
% figure;stem( HF_MI_info_fL_final(:,2) , HF_MI_info_fL_final(:,5) );
%
%
% figure;plot3( IAmp_ana_dat(:,2) , IAmp_ana_dat(:,3) , IAmp_ana_dat(:,6) , 'o');

%{
fL_ind_lgc = HF_MI_info(:,1) == 9;
HF_MI_info_fL = HF_MI_info(fL_ind_lgc, [1,2,5]); HF_MI_info_fL = HF_MI_info_fL( HF_MI_info_fL(:,3) ~= 0 , : );
figure; stem( HF_MI_info_fL(:,2) , HF_MI_info_fL(:,3) );
ip = findchangepts(HF_MI_info_fL(:,3), "Statistic","std" , MinThreshold=2);
HF_MI_info_fL(ip,2)
%}
end

%%




function [ipt_freq, HF_MI_info_fL_1, HF_MI_info_fL_final, break_loop] = segmt_mean(HF_MI_info, fL_ind_lgc, no_fH_vals)
break_loop = 0;

%Non-zeros items are used to calculate the averages for segmentation
HF_MI_info_fL = HF_MI_info(fL_ind_lgc, [1,2,5]); %HF_MI_info_fL = HF_MI_info_fL( HF_MI_info_fL(:,3) ~= 0 , : );

HF_MI_info_fL_final = HF_MI_info(fL_ind_lgc, [1,2,5]);
HF_MI_info_fL_1 = HF_MI_info_fL_final;

ipt = findchangepts(HF_MI_info_fL(:,3), "Statistic", "linear", "MinThreshold", 0.5);
ipt_ind = [1;ipt;length(HF_MI_info_fL)]; %HF_MI_info_fL(ipt_ind, [2,3])
if isempty(ipt) || length(ipt)<=2
    break_loop = 1; 
    HF_MI_info_fL = 0;
    ipt_freq = 0;
    return;
else
    ipt_freq = HF_MI_info_fL(ipt_ind , 2);
end

%{
min_max = [min(HF_MI_info_fL(:,3)) , max(HF_MI_info_fL(:,3))];
figure;plot( HF_MI_info_fL(:,2) , HF_MI_info_fL(:,3) ); hold on; grid on;
for i=1:length(ipt); plot( HF_MI_info_fL(ipt(i),2).*ones(1,10+1) , [min_max(1):(min_max(2)-min_max(1))/10:min_max(2)] , 'r' , LineWidth=1);end
sgtitle( num2str(HF_MI_info_fL(1,1)) );
%}

HF_MI_seg_mean = zeros( no_fH_vals , 1 );
for i=1:length(ipt_ind)-1
    ind_rng = [ find( HF_MI_info_fL_final(:,2) == HF_MI_info_fL(ipt_ind(i) , 2) ) , find( HF_MI_info_fL_final(:,2) == HF_MI_info_fL(ipt_ind(i+1) , 2) ) ];
    ind_rng = ind_rng(1):ind_rng(2);
    HF_MI_seg_mean(ind_rng,1) = mean( HF_MI_info_fL( ipt_ind(i):ipt_ind(i+1) , 3) ); %figure;stem( HF_MI_info_fL(:,2) , HF_MI_seg_mean );
end
HF_MI_seg_mean( HF_MI_seg_mean<0 & HF_MI_seg_mean>-0.1 ,:) = 0;
HF_MI_seg_mean( HF_MI_seg_mean>0  & HF_MI_seg_mean<0.1 ,:) = 0;

HF_MI_seg_mean( HF_MI_seg_mean<0 ,:) = -1;
HF_MI_seg_mean( HF_MI_seg_mean>0 ,:) = 1;

HF_MI_info_fL_final = HF_MI_info(fL_ind_lgc, :);
HF_MI_info_fL_final(:,end) = HF_MI_seg_mean;

HF_MI_info_fL_1 = HF_MI_info_fL_final(:, [1,2,5]);

HF_MI_info_fL_final = [HF_MI_info_fL_final , zeros(no_fH_vals,1)];
for i=1:no_fH_vals
    if HF_MI_info_fL_final(i,3) ~= 0
        pos_lgc = HF_MI_info_fL_final(:,2) == HF_MI_info_fL_final(i,3);
        if sum(pos_lgc) ~= 0
            HF_MI_info_fL_final(i,6) = HF_MI_info_fL_final(pos_lgc,5);
        end
    else
        break_loop = 1; return;
    end
end
HF_MI_info_fL_final = [HF_MI_info_fL_final , zeros(no_fH_vals,1)];
for i=1:no_fH_vals
    if HF_MI_info_fL_final(i,4) ~= 0
        pos_lgc = HF_MI_info_fL_final(:,2) == HF_MI_info_fL_final(i,4);
        if sum(pos_lgc) ~= 0
            HF_MI_info_fL_final(i,7) = HF_MI_info_fL_final(pos_lgc,5);
        end
    else
        break_loop = 1; return;
    end
end

end

%%

function [HF_MI_info_fL_1, HF_MI_info_fL_final, break_loop] = segmt_mean_2(HF_MI_info, fL_ind_lgc, no_fH_vals)
break_loop = 0;
HF_MI_info_fL_final = HF_MI_info(fL_ind_lgc, :);

HF_MI_info_fL_1 = HF_MI_info_fL_final(:, [1,2,5]);

HF_MI_info_fL_final = [HF_MI_info_fL_final , zeros(no_fH_vals,1)];
for i=1:no_fH_vals
    if HF_MI_info_fL_final(i,3) ~= 0
        pos_lgc = HF_MI_info_fL_final(:,2) == HF_MI_info_fL_final(i,3);
        if sum(pos_lgc) ~= 0
            HF_MI_info_fL_final(i,6) = HF_MI_info_fL_final(pos_lgc,5);
        end
    else
        break_loop = 1; return;
    end
end
HF_MI_info_fL_final = [HF_MI_info_fL_final , zeros(no_fH_vals,1)];
for i=1:no_fH_vals
    if HF_MI_info_fL_final(i,4) ~= 0
        pos_lgc = HF_MI_info_fL_final(:,2) == HF_MI_info_fL_final(i,4);
        if sum(pos_lgc) ~= 0
            HF_MI_info_fL_final(i,7) = HF_MI_info_fL_final(pos_lgc,5);
        end
    else
        break_loop = 1; return;
    end
end

end