function [Comod_harmonic] = AM_var_commod(All_freq_comb_1, phs_data_mat, fL_vals, fH_vals)
%AM_VAR_COMMOD Quantify amplitude-modulation variation for detected PAC pairs.
%   This diagnostic uses the noise-free canonical model components stored in
%   PHS_DATA_MAT. For each detected low/high-frequency pair, it computes the
%   standard-deviation ratio between the summed high-frequency contribution
%   and the low-frequency contribution, then maps that value to the requested
%   frequency grid.
%   This routine supports post-identification inspection of the model
%   decomposition described in Sections III C and III F of the paper.
%   Inputs
%   ------
%   All_freq_comb_1 : Numeric matrix whose first two columns contain the
%                     detected low- and high-frequency centres in Hz.
%   phs_data_mat     : Cell array aligned with All_freq_comb_1. Each cell
%                     contains simulated canonical components; column 1 is
%                     the low-frequency response and columns 2:end are the
%                     high-frequency linear and interaction responses.
%   fL_vals          : Vector of queried low-frequency centres in Hz.
%   fH_vals          : Vector of queried high-frequency centres in Hz.
%   Output
%   ------
%   Comod_harmonic   : numel(fH_vals)-by-numel(fL_vals) map containing the
%                     amplitude-variation ratio at detected pairs and zero
%                     elsewhere.

%% Compute the diagnostic value for each detected frequency pair

pos_freq_comp_vec = All_freq_comb_1(:,[1,2]);
no_freqs = size(All_freq_comb_1,1);
IF_harmonic_test_dat = zeros(no_freqs , 3);

for i = 1:no_freqs
    probe_freq = All_freq_comb_1(i ,[1,2]);
    loc_prb_freq = sum(pos_freq_comp_vec == probe_freq,2)==2;
    Phs_dat_prb_freq = phs_data_mat{loc_prb_freq, 1}; 
    IF_harmonic_test_dat(i,:) = [ probe_freq , ( std( sum( Phs_dat_prb_freq(:,2:end) ) ) / std( Phs_dat_prb_freq(:,1) ) ) ];
end

%% Map the pairwise values to the full comodulogram grid
[fL_grd, fH_grd] = meshgrid(fL_vals, fH_vals);
pos_freq_comp_vec = [fL_grd(:) , fH_grd(:)];
no_probes = size(pos_freq_comp_vec,1);

harmonic_freq = zeros(no_probes,1);
for i = 1:no_freqs
    harmonic_freq(sum(pos_freq_comp_vec == IF_harmonic_test_dat(i ,[1,2]), 2) == 2, 1) =  IF_harmonic_test_dat(i , 3) ;
end
Comod_harmonic = reshape( harmonic_freq , size(fL_grd));


end
