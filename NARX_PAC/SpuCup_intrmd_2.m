function [HF_MI_score_1, HF_MI_score_final, Comod_intrmd_rmv] = SpuCup_intrmd_2(fL_vals, fH_vals, Comod, diff_comod, All_freq_comb_1, phs_data_mat, Fs, sine_freq)
%SPUCUP_INTRMD_2 Suppress intermodulation-related spurious PAC detections.
%   The routine evaluates the local alternating-sign structure of the
%   discriminator D = sign(H_norm - C_norm), where H is the simulated
%   high-frequency magnitude map and C is the NARX-PAC modulation-index map.
%   A candidate high-frequency component is supported when its discriminator
%   sign is opposite to the signs at fH-fL and fH+fL. This follows the
%   discriminator construction in Section III E, equations (21)-(23).
%   For broadband/non-stationary data (sine_freq == 0), D is first segmented
%   along the high-frequency axis and each segment is represented by the sign
%   of its mean. For sinusoidal test data, the unsegmented binwise signs are
%   used directly.
%   Inputs
%   ------
%   fL_vals, fH_vals  : Low- and high-frequency query vectors in Hz.
%   Comod              : NARX-PAC modulation-index comodulogram.
%   diff_comod         : Unthresholded discriminator precursor H_norm-C_norm.
%   All_freq_comb_1    : Retained interface argument containing detected
%                        frequency-pair information; not used by the current
%                        scoring implementation.
%   phs_data_mat        : Retained interface argument containing canonical
%                        simulations; only its dimensions are read currently.
%   Fs                  : Retained sampling-frequency argument in Hz.
%   sine_freq           : Processing selector. Use 0 for segmented broadband
%                        processing and a non-zero value for direct binwise
%                        processing.
%   Outputs
%   -------
%   HF_MI_score_1       : [fL, fH, score] for every grid point. Negative scores
%                        mark intermodulation-related detections for removal.
%   HF_MI_score_final   : Frequency-pair table [fL, fH, fH-fL, fH+fL, 0]
%                        retained for compatibility with existing workflows.
%   Comod_intrmd_rmv    : Comod after suppressing grid points with a negative
%                        intermodulation-discrimination score.

%% Assemble the discriminator table for the frequency grid

%% MI and HF data

no_fH_vals = length(fH_vals);
N = size( phs_data_mat{1, 1} , 1);
% Construct the complete low/high-frequency grid.
[fL_grd, fH_grd] = meshgrid(fL_vals, fH_vals);
pos_freq_comp_vec = [fL_grd(:) , fH_grd(:)];

HF_MI_info = [ pos_freq_comp_vec , ( pos_freq_comp_vec(:,2) + [-pos_freq_comp_vec(:,1) , pos_freq_comp_vec(:,1) ] ) ];


comod_D = diff_comod;
comod_D( comod_D > 0 ) = 1; comod_D( comod_D < 0 ) = -1; 

HF_MI_info = [ HF_MI_info , comod_D(:) ]; 


% Preserve the complete grid description for downstream/legacy consumers.
HF_MI_score_final = HF_MI_info; HF_MI_score_final(:,end) = 0;
%% MI and HF info per low frequency
HF_MI_score_1 = zeros( size(HF_MI_score_final,1) , 3 );
HF_MI_info_2 = cell( length(fL_vals) , 1 );
commod_intrmd_HF_clstrs = cell( length(fL_vals) , 1 );

%% Process each low-frequency slice independently
loop_cnt = 1;
for fL_freq = fL_vals
    fL_ind_lgc = HF_MI_info(:,1) == fL_freq;

    break_loop = 0;
    %% Segmentation and scoring -- MIHF

    % Segmentation reduces isolated sign fluctuations in broadband cases.
    if sine_freq == 0

        [ipt_freq, HF_MI_info_fL_1, HF_MI_info_fL_final, break_loop] = segmt_mean(HF_MI_info, fL_ind_lgc, no_fH_vals);


        HF_MI_info_2{loop_cnt,1} = HF_MI_info_fL_1;
        commod_intrmd_HF_clstrs{loop_cnt,1} = ipt_freq;

    else
        [~, HF_MI_info_fL_final, break_loop] = segmt_mean_2(HF_MI_info, fL_ind_lgc, no_fH_vals);
    end
    if break_loop == 1
        HF_MI_score_1(fL_ind_lgc,:) = [ HF_MI_info_fL_final(:, [1,2]) , 0.*HF_MI_info_fL_final(:,end)];
        loop_cnt = loop_cnt + 1;
        continue;
    end
    % Positive evidence
    % Positive evidence is obtained when D at the candidate component has
    % the opposite sign to D at either immediate intermodulation frequency.
    pos_evd_Idiff = max( [ zeros(no_fH_vals,1), -sign( HF_MI_info_fL_final(:,5) ) .* HF_MI_info_fL_final(:,6) ] , [] , 2 );
    pos_evd_Ism = max( [ zeros(no_fH_vals,1), -sign( HF_MI_info_fL_final(:,5) ) .* HF_MI_info_fL_final(:,7) ] , [] , 2  );
    pos_evd = pos_evd_Idiff + pos_evd_Ism;

    % Absolute Negative evidence
    % Frequencies identified as immediate intermodulations of a supported
    % candidate receive explicit negative evidence.
    neg_evd_HFs = sort( unique( [ HF_MI_info_fL_final( pos_evd == 2 , 3) ; HF_MI_info_fL_final( pos_evd == 2 , 4) ] ) , 'ascend');
    neg_evd = zeros(no_fH_vals,1);
    for i=1:length(neg_evd_HFs); neg_evd( HF_MI_info_fL_final(:,2) == neg_evd_HFs(i) , 1) = -1; end
    
    % Negative evidence--Intermodulations in the edges
    neg_evd_edge = -( pos_evd_Idiff.*(HF_MI_info_fL_final(:,7)==0) + pos_evd_Ism.*( HF_MI_info_fL_final(:,6)==0) );
    
    pos_evd_edge_Ism_HFs = HF_MI_info_fL_final( (pos_evd_Idiff.*(HF_MI_info_fL_final(:,7)==0))==1 , 3);
    pos_evd_edge_Ism = zeros(no_fH_vals,1);
    for i=1:length(pos_evd_edge_Ism_HFs); pos_evd_edge_Ism( HF_MI_info_fL_final(:,2) == pos_evd_edge_Ism_HFs(i) , 1) = 1; end
    pos_evd_edge_Idiff_HFs = HF_MI_info_fL_final( (pos_evd_Ism.*(HF_MI_info_fL_final(:,6)==0))==1 , 4);
    pos_evd_edge_Idiff = zeros(no_fH_vals,1);
    for i=1:length(pos_evd_edge_Idiff_HFs); pos_evd_edge_Idiff( HF_MI_info_fL_final(:,2) == pos_evd_edge_Idiff_HFs(i) , 1) = 1; end
    pos_evd_edge = pos_evd_edge_Idiff + pos_evd_edge_Ism;

    score = 0.5.*pos_evd + 0.*pos_evd_edge + neg_evd + 0.*neg_evd_edge;


    HF_MI_score_1(fL_ind_lgc,:) = [ HF_MI_info_fL_final(:, [1,2]) , score];

    
loop_cnt = loop_cnt + 1;
end

%% Convert pairwise scores to a mask and apply it to the comodulogram
Comod_intrmd = reshape( HF_MI_score_1(:,3) , size(fL_grd)) >= 0;
Comod_intrmd_rmv = Comod .* Comod_intrmd;


end


function [ipt_freq, HF_MI_info_fL_1, HF_MI_info_fL_final, break_loop] = segmt_mean(HF_MI_info, fL_ind_lgc, no_fH_vals)
%SEGMT_MEAN Segment one low-frequency discriminator slice.
%   Change points are estimated along the high-frequency axis. Each segment
%   is replaced by the sign of its mean discriminator value, after which the
%   discriminator signs at fH-fL and fH+fL are appended for scoring.
%   Inputs: HF_MI_info is [fL, fH, fH-fL, fH+fL, D]; fL_ind_lgc selects one
%   low-frequency slice; no_fH_vals is the number of high-frequency bins.
%   Outputs: ipt_freq contains segment boundaries in Hz; HF_MI_info_fL_1 is
%   a compact [fL, fH, segmentedD] table; HF_MI_info_fL_final appends the D
%   values at both immediate intermodulations; break_loop flags slices that
%   cannot be evaluated reliably.
break_loop = 0;

% Use the discriminator trace to estimate broad sign-consistent segments.
HF_MI_info_fL = HF_MI_info(fL_ind_lgc, [1,2,5]); %HF_MI_info_fL = HF_MI_info_fL( HF_MI_info_fL(:,3) ~= 0 , : );

HF_MI_info_fL_final = HF_MI_info(fL_ind_lgc, [1,2,5]);
HF_MI_info_fL_1 = HF_MI_info_fL_final;

% Detect broad sign regions rather than scoring individual noisy bins.
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

% Append the discriminator value at the lower intermodulation fH-fL.
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
% Append the discriminator value at the upper intermodulation fH+fL.
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


function [HF_MI_info_fL_1, HF_MI_info_fL_final, break_loop] = segmt_mean_2(HF_MI_info, fL_ind_lgc, no_fH_vals)
%SEGMT_MEAN_2 Prepare an unsegmented discriminator slice.
%   This variant keeps the original binwise discriminator signs and appends
%   the signs at fH-fL and fH+fL. It is used for sinusoidal test cases where
%   broad frequency segmentation is unnecessary.
%   Inputs and outputs follow SEGMT_MEAN, except that no change-point
%   frequencies are returned.
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
