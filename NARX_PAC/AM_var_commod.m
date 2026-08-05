function [Comod_harmonic] = AM_var_commod(All_freq_comb_1, phs_data_mat, fL_vals, fH_vals)

pos_freq_comp_vec = All_freq_comb_1(:,[1,2]);
no_freqs = size(All_freq_comb_1,1);
IF_harmonic_test_dat = zeros(no_freqs , 3);

for i = 1:no_freqs
    probe_freq = All_freq_comb_1(i ,[1,2]);
    loc_prb_freq = sum(pos_freq_comp_vec == probe_freq,2)==2;
    Phs_dat_prb_freq = phs_data_mat{loc_prb_freq, 1}; 
    IF_harmonic_test_dat(i,:) = [ probe_freq , ( std( sum( Phs_dat_prb_freq(:,2:end) ) ) / std( Phs_dat_prb_freq(:,1) ) ) ];
end

[fL_grd, fH_grd] = meshgrid(fL_vals, fH_vals);
pos_freq_comp_vec = [fL_grd(:) , fH_grd(:)];
no_probes = size(pos_freq_comp_vec,1);

harmonic_freq = zeros(no_probes,1);
for i = 1:no_freqs
    harmonic_freq(sum(pos_freq_comp_vec == IF_harmonic_test_dat(i ,[1,2]), 2) == 2, 1) =  IF_harmonic_test_dat(i , 3) ;
end
Comod_harmonic = reshape( harmonic_freq , size(fL_grd));
%-------------


% figure; imagesc(fL_vals, fH_vals, Comod_harmonic); colorbar; axis xy; set(gca, 'FontSize', 18);

end