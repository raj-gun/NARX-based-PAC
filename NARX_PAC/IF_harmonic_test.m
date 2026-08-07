function [Comod_harmonic_rmv, IF_harmonic_test_dat] = IF_harmonic_test (All_freq_comb_1, phs_data_mat, fL_vals, fH_vals, Comod, Ts)
%IF_HARMONIC_TEST Reject harmonic-related spurious PAC detections.
%   For each detected pair, the canonical model components are summed and
%   converted to an analytic signal. The instantaneous frequency is obtained
%   from the time derivative of the unwrapped analytic phase. For numerical
%   implementation, a 95% non-negative instantaneous-frequency threshold is
%   used to tolerate small deviations introduced by numerical phase estimation
%   and differentiation, while approximating the theoretical condition in
%   Section III D, equation (18), of the paper.
%   Inputs
%   ------
%   All_freq_comb_1      : Numeric matrix whose first two columns contain the
%                          detected low- and high-frequency centres in Hz.
%   phs_data_mat          : Cell array of noise-free canonical component
%                          simulations aligned with All_freq_comb_1.
%   fL_vals, fH_vals      : Low- and high-frequency query vectors in Hz.
%   Comod                  : NARX-PAC modulation-index comodulogram.
%   Ts                     : Sampling interval in seconds.
%   Outputs
%   -------
%   Comod_harmonic_rmv    : Comod with harmonic-related detections set to zero.
%   IF_harmonic_test_dat  : One row per detected pair: [fL, fH, harmonicFlag],
%                          where harmonicFlag is 1 for a rejected pair.

%% Evaluate the instantaneous-frequency criterion for detected pairs

pos_freq_comp_vec = All_freq_comb_1(:,[1,2]);
no_freqs = size(All_freq_comb_1,1);
IF_harmonic_test_dat = zeros(no_freqs , 3);

for i = 1:no_freqs
    probe_freq = All_freq_comb_1(i ,[1,2]);
    loc_prb_freq = sum(pos_freq_comp_vec == probe_freq,2)==2;
    Phs_dat_prb_freq = phs_data_mat{loc_prb_freq, 1}; 
    IFreq = diff( unwrap( angle( hilbert( sum(Phs_dat_prb_freq,2) ) ) ) ) ./ Ts;
    % Use a 95% threshold as a numerical tolerance for small local deviations
    % introduced by discrete Hilbert-phase estimation and differentiation; the
    % theoretical criterion is positive instantaneous frequency for all t.
    IF_cond = IFreq >= 0;
    IF_harmonic_test_dat(i,:) = [ probe_freq , ( sum( IF_cond ) / (size(Phs_dat_prb_freq,1)-1) ) >=0.95 ];
end

%% Construct the rejection mask on the complete frequency grid
[fL_grd, fH_grd] = meshgrid(fL_vals, fH_vals);
pos_freq_comp_vec = [fL_grd(:) , fH_grd(:)];
no_probes = size(pos_freq_comp_vec,1);

harmonic_freq = zeros(no_probes,1);
for i = 1:no_freqs
    harmonic_freq(sum(pos_freq_comp_vec == IF_harmonic_test_dat(i ,[1,2]), 2) == 2, 1) =  IF_harmonic_test_dat(i , 3) ;
end
Comod_harmonic = reshape( harmonic_freq , size(fL_grd)) ~= 1;
Comod_harmonic_rmv = Comod .* Comod_harmonic;


end
