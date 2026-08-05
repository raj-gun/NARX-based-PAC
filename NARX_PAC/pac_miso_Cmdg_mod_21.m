function [Comods , diff_comod, phs_data_mat_2, fL_grd, fH_grd, All_freq_comb_2, mod_trm_clstr_ERR_2, All_freq_comb_ARX_1, All_freq_comb_ARX_2, narx_pac_modls_1 , narx_pac_modls_2] = pac_miso_Cmdg_mod_21(signal, fL_vals, fH_vals, Fs, RCT, filt_typ, frq_bndw_LF, frq_bndw_HF)

Ts = 1/Fs;
N = length(signal);
fft_res = Fs/N;

% Generate the grid
[fL_grd, fH_grd] = meshgrid(fL_vals, fH_vals);
pos_freq_comp_vec = [fL_grd(:) , fH_grd(:)];
no_probes = size(pos_freq_comp_vec,1);

std_cos = sqrt(0.5); signal = ( (signal-mean(signal)).* ( std_cos/std(signal) ) );

phi = [0 , 0];
vrb = 0;

coupling_coeff_arx = zeros(no_probes,5);
% arx_pac_modls = cell(no_probes,1);

coupling_coeff = zeros(no_probes,5);
narx_pac_modls = cell(no_probes,1);
phs_data_mat = cell(no_probes,1);

intrmd_ratio_thrshld = 0.7;
HF_LF_cond = [0,0.7];
intrmd_high_freq_diff_thrshld = 0.2;
MI_cond = [0.2, 5];

%% Extract narrow-band frequency data

K_LF = round(frq_bndw_LF/fft_res); if K_LF==0; K_LF=1; end %floor( ( ( (frq_bndw_LF)/Fs )*N ) );
K_HF = round(frq_bndw_HF/fft_res); if K_HF==0; K_HF=1; end %floor( ( ( (frq_bndw_HF)/Fs )*N ) );

[fL_bndpss_dat] = freqs_filt_dat(fL_vals, K_LF, signal, Fs, N, filt_typ{1});
[fH_bndpss_dat] = freqs_filt_dat(fH_vals, K_HF, signal, Fs, N, filt_typ{2});

[fL_grd_ind, fH_grd_ind] = meshgrid( 1:length(fL_vals), 1:length(fH_vals) );
pos_freq_comp_ind = [fL_grd_ind(:) , fH_grd_ind(:)];
freqs_filt_dat_cntnr = zeros(N,2,no_probes);%cell(no_probes,1);
for i=1:no_probes
    fL_ind = pos_freq_comp_ind(i,1);
    fH_ind = pos_freq_comp_ind(i,2);
    freqs_filt_dat_cntnr(:,:,i) = [fL_bndpss_dat(:,fL_ind) , fH_bndpss_dat(:,fH_ind)];
end
%% Scan freq grid with ARX
% ==================================================================
% Scan the grid
% ==================================================================
parfor i = 1:no_probes % find(sum(pos_freq_comp_vec==[7,63],2)==2) %
    pos_freq_comp = pos_freq_comp_vec(i,:);
    s_LF_HF_filt_dat = freqs_filt_dat_cntnr(:,:,i); % freqs_filt_dat_cntnr{i,1}; 2

    [cplng_param,~] = grd_srch_lin(signal, s_LF_HF_filt_dat, Fs, Ts, pos_freq_comp, phi, vrb, RCT);
    coupling_coeff_arx(i,:) = cplng_param;
end
% ==================================================================
% ==================================================================

% ==================================================================
coupling_coeff_arx(isnan(coupling_coeff_arx)) = 0;
coupling_coeff_arx( repmat( coupling_coeff_arx(:,1)==0 | coupling_coeff_arx(:,2)==0 , [1,size(coupling_coeff_arx,2)] ) ) = 0;

logic_ind_cplcff = coupling_coeff_arx(:,1)~=0 & coupling_coeff_arx(:,2)~=0;

All_freq_comb_ARX_1 = [pos_freq_comp_vec(logic_ind_cplcff , :) , coupling_coeff_arx(logic_ind_cplcff , :)];
intrmd_low_freq_ratio = All_freq_comb_ARX_1(:,[3,4]) ./ All_freq_comb_ARX_1(:,[6,6]);
intrmd_low_high_freq_ratio = All_freq_comb_ARX_1(:,[3,4]) ./ ( All_freq_comb_ARX_1(:,[6,6]) .* All_freq_comb_ARX_1(:,[5,5]) );
All_freq_comb_ARX_1 = [ All_freq_comb_ARX_1 , intrmd_low_freq_ratio , intrmd_low_high_freq_ratio];

HF_LF_cond_lgc_1 = ( All_freq_comb_ARX_1(:,5)./All_freq_comb_ARX_1(:,6) ) <= HF_LF_cond(2);
HF_LF_cond_lgc_2 = ( All_freq_comb_ARX_1(:,5)./All_freq_comb_ARX_1(:,6) ) >= HF_LF_cond(1);
HF_LF_cond_lgc = HF_LF_cond_lgc_1 & HF_LF_cond_lgc_2;

pos_freq_comp_vec_ind = 1:no_probes;
All_freq_comb_ARX_pos_freq_ind_incld = pos_freq_comp_vec_ind(logic_ind_cplcff);
Final_pos_freq_ind_incld = All_freq_comb_ARX_pos_freq_ind_incld;%(HF_LF_cond_lgc);
% ==================================================================


%% Identify PAC frequency pairs
no_probes_1 = length(Final_pos_freq_ind_incld);
coupling_coeff_1 = zeros(no_probes_1,5);
phs_data_mat_1 = cell(no_probes_1,1); %cell(1,no_probes_1);%
narx_pac_modls_1 = cell(no_probes_1,1);

pos_freq_comp_vec_1 = pos_freq_comp_vec(Final_pos_freq_ind_incld,:);
freqs_filt_dat_cntnr_1 = freqs_filt_dat_cntnr(:,:,Final_pos_freq_ind_incld);

% ==================================================================
% Identify PAC frequency pairs
% ==================================================================
% for i = find(sum(pos_freq_comp_vec_1==[7,57],2)==2) %
parfor i = 1:no_probes_1 % find(sum(pos_freq_comp_vec_1==[7,63],2)==2) %
    %try
    pos_freq_comp = pos_freq_comp_vec_1(i,:);
    s_LF_HF_filt_dat = freqs_filt_dat_cntnr_1(:,:,i);%freqs_filt_dat_cntnr{i,1};2

    [cplng_param, model, phs_data] = grd_srch(signal, s_LF_HF_filt_dat, Fs, Ts, pos_freq_comp, phi, vrb, RCT);
    coupling_coeff_1(i,:) = cplng_param;
    phs_data_mat_1{i,1} = phs_data;
    narx_pac_modls_1{i,1} = model;
    %{
    catch
    disp(['innerloop=',  num2str(i)]);
    end
    %}
end
coupling_coeff(Final_pos_freq_ind_incld,:) = coupling_coeff_1;
for j = 1:no_probes_1
    narx_pac_modls{Final_pos_freq_ind_incld(j),1} = narx_pac_modls_1{j,1};
    phs_data_mat{ Final_pos_freq_ind_incld(j) , :} = phs_data_mat_1{j,:};
end
% ==================================================================
% ==================================================================

%% ==================================================================
%% Process MISO NARX results
%% ==================================================================

%-----------------------------------------------------------------------
coupling_coeff(isnan(coupling_coeff)) = 0;
coupling_coeff(isinf(coupling_coeff)) = 0;
coupling_coeff( repmat( coupling_coeff(:,1)==0 | coupling_coeff(:,2)==0 , [1,size(coupling_coeff,2)] ) ) = 0;

logic_ind_cplcff = coupling_coeff(:,1)~=0 & coupling_coeff(:,2)~=0;

if sum(logic_ind_cplcff) ~= 0
    narx_pac_modls_1 = narx_pac_modls(logic_ind_cplcff , 1);
    phs_data_mat_1 = phs_data_mat(logic_ind_cplcff, 1);
    %
    mod_trm_clstr_ERR_1 = zeros(length(narx_pac_modls_1) , 6);
    for i = 1:length(narx_pac_modls_1)
        match_ind = intr_modul_comp(narx_pac_modls_1{i,1});
        ERR_u1u2 = sum( narx_pac_modls_1{i,1}{1,16}.ERR(logical(match_ind)) );
        %
        match_ind = lin_comp_u2(narx_pac_modls_1{i,1});
        ERR_u2 = sum( narx_pac_modls_1{i,1}{1,16}.ERR(logical(match_ind)) );
        %
        match_ind = lin_comp_u1(narx_pac_modls_1{i,1});
        ERR_u1 = sum( narx_pac_modls_1{i,1}{1,16}.ERR(logical(match_ind)) );
        %
        mod_trm_clstr_ERR_1(i,:) = [ERR_u1 , ERR_u2 , ERR_u1u2 , ERR_u1u2/ERR_u2 , ERR_u2/ERR_u1 , ERR_u1u2/(ERR_u1+ERR_u2)];
    end
    mod_trm_clstr_ERR_1 = [pos_freq_comp_vec(logic_ind_cplcff , :) , mod_trm_clstr_ERR_1];
else
    narx_pac_modls_1 = 'N\A';
    mod_trm_clstr_ERR_1 = 'N\A';
end

All_freq_comb_1 = [pos_freq_comp_vec(logic_ind_cplcff , :) , coupling_coeff(logic_ind_cplcff , :)];
intrmd_low_freq_ratio = All_freq_comb_1(:,[3,4]) ./ All_freq_comb_1(:,[6,6]);
intrmd_low_high_freq_ratio = All_freq_comb_1(:,[3,4]) ./ ( All_freq_comb_1(:,[6,6]) .* All_freq_comb_1(:,[5,5]) );
All_freq_comb_1 = [ All_freq_comb_1 , intrmd_low_freq_ratio , intrmd_low_high_freq_ratio];
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
MI_cond_lgc = ( All_freq_comb_1(:,7) > MI_cond(1) ) & ( All_freq_comb_1(:,7) < MI_cond(2) );
%-------------
HF_LF_cond_lgc_1 = ( All_freq_comb_1(:,5)./All_freq_comb_1(:,6) ) <= HF_LF_cond(2);
HF_LF_cond_lgc_2 = ( All_freq_comb_1(:,5)./All_freq_comb_1(:,6) ) >= HF_LF_cond(1);
HF_LF_cond_lgc = HF_LF_cond_lgc_1 & HF_LF_cond_lgc_2;
%-------------
% intrmd_freqs = All_freq_comb_1(:,[3,4]);
% intrmd_cond_lgc = ( min(intrmd_freqs,[],2) ./ max(intrmd_freqs,[],2) ) > intrmd_ratio_thrshld;
%-------------
% modltn_cond = All_freq_comb_1(:,7) < 0.8 | All_freq_comb_1(:,7) > 1.2;
%-------------
% intrmd_high_freq_diff = ...
% [abs( All_freq_comb_1(:,5) - All_freq_comb_1(:,3) ) ./ max(All_freq_comb_1(:,5) , All_freq_comb_1(:,3)) , ...
% abs( All_freq_comb_1(:,5) - All_freq_comb_1(:,4) ) ./ max(All_freq_comb_1(:,5) , All_freq_comb_1(:,4))] > intrmd_high_freq_diff_thrshld;
% intrmd_high_freq_diff = sum(intrmd_high_freq_diff,2) == 2;
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
% overall_logic = intrmd_cond_lgc & intrmd_high_freq_diff;
% overall_logic = HF_LF_cond_lgc & intrmd_cond_lgc;
% overall_logic = (HF_LF_cond_lgc & intrmd_cond_lgc) & intrmd_high_freq_diff;
overall_logic = HF_LF_cond_lgc & MI_cond_lgc;
%-------------
All_freq_comb_2 = All_freq_comb_1( overall_logic, :);
narx_pac_modls_2 = narx_pac_modls_1( overall_logic, : );
mod_trm_clstr_ERR_2 = mod_trm_clstr_ERR_1( overall_logic, : );
% spect_comp_mat_2 = spect_comp_mat_1(overall_logic,:);
All_freq_comb_ARX_2 = All_freq_comb_ARX_1( overall_logic, : );
phs_data_mat_2 = phs_data_mat_1( overall_logic, : );
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
coupling_coeff_raw = zeros(no_probes,1);
for i = 1:size(All_freq_comb_1,1)
    coupling_coeff_raw(sum(pos_freq_comp_vec == All_freq_comb_1(i ,[1,2]), 2) == 2, 1) =  All_freq_comb_1(i ,7) ;
end
Comod_raw = reshape(coupling_coeff_raw,size(fL_grd));
%-------------
HF_vals = zeros(no_probes,1);
for i = 1:size(mod_trm_clstr_ERR_2,1)
    HF_vals(sum(pos_freq_comp_vec == mod_trm_clstr_ERR_2(i ,[1,2]), 2) == 2, 1) =  All_freq_comb_2(i , 5); % mod_trm_clstr_ERR_2(i , 4); %
end
Comod_HF = reshape(HF_vals,size(fL_grd));
%-------------
MI_ERR_vals = zeros(no_probes,1);
for i = 1:size(mod_trm_clstr_ERR_2,1)
    MI_ERR_vals(sum(pos_freq_comp_vec == mod_trm_clstr_ERR_2(i ,[1,2]), 2) == 2, 1) =  mod_trm_clstr_ERR_2(i , 6);
end
Comod_MI_ERR = reshape(MI_ERR_vals,size(fL_grd));
%-------------
coupling_coeff_refine = zeros(no_probes,1);
for i = 1:size(All_freq_comb_2,1)
    coupling_coeff_refine(sum(pos_freq_comp_vec == All_freq_comb_2(i ,[1,2]), 2) == 2, 1) =  All_freq_comb_2(i ,7);
end
Comod = reshape(coupling_coeff_refine,size(fL_grd));
%-----------------------------------------------------------------------

Comod_norm = Comod./std(Comod,0,1); Comod_norm(isnan(Comod_norm)) = 0;
Comod_HF_norm = Comod_HF./std(Comod_HF,0,1); Comod_HF_norm(isnan(Comod_HF_norm)) = 0;
diff_comod = Comod_HF_norm - Comod_norm;

Comods = {Comod,Comod_raw};
%%

%% ==================================================================
%% MISO NARX surrogate results
%% ==================================================================

end

%% Local functions

%-------------------------------
function [freq_bndpss_dat] = freqs_filt_dat(freq_vals, frq_bndw, signal, Fs, N, filt_typ)
n_freq_vals = length(freq_vals);
freq_bndpss_dat = zeros(N, n_freq_vals);
for i=1:length(freq_vals)
    sig_filt = nrrw_bnd_fft_filt(signal', Fs, freq_vals(i), frq_bndw, filt_typ);
    if size(sig_filt,2) ~= 1; sig_filt = sig_filt.'; end % Transpose to a column vector
    freq_bndpss_dat(:,i) = sig_filt;
end
end
%-------------------------------

%-------------------------------
function [cplng_param, model, phs_data] = grd_srch(signal, env, Fs, Ts, pos_freq_comp, phi, vrb, RCT)
freq_pos = @(freq,Fs,fftn) floor(freq*fftn/Fs)+1;
% freq_pos = @(freq,Fs,fftn) round(freq*fftn/Fs);
if pos_freq_comp(2) > pos_freq_comp(1)
    nb2 = [ceil( (1/pos_freq_comp(1))/(4*Ts) ) , ceil( (1/pos_freq_comp(2))/(Ts) )];
    [mods_nofrfs, model] = pac_miso_base(nb2, signal, env, Fs, Ts, pos_freq_comp, phi, vrb, RCT);
    fft_mags = mods_nofrfs{1,9};
    fftn = length(fft_mags);
    phs_data = mods_nofrfs{1,12};

    intrmd_sm = pos_freq_comp(2)+pos_freq_comp(1); intrmd_dff = pos_freq_comp(2)-pos_freq_comp(1);

    low_freq_mag  = fft_mags(freq_pos(pos_freq_comp(1),Fs,fftn),:);
    high_freq_mag = fft_mags(freq_pos(pos_freq_comp(2),Fs,fftn),:);
    intrmd_sm_freq_mag  = fft_mags(freq_pos(intrmd_sm,Fs,fftn),:);
    intrmd_dff_freq_mag  = fft_mags(freq_pos(intrmd_dff,Fs,fftn),:);

	intrmd_freqs = [intrmd_sm_freq_mag,intrmd_dff_freq_mag];
	intrmd_high_freq_ratio = intrmd_freqs ./ high_freq_mag;
	modl = mean( intrmd_high_freq_ratio );
	cplng_param = [ intrmd_dff_freq_mag , intrmd_sm_freq_mag , high_freq_mag , low_freq_mag , modl ];
else
    cplng_param = [0,0,0,0,0];
    model = 0;
    phs_data = 0;
end
% fftn_1 = length(signal);
% signal_fft = (Ts/length(signal)).*fft(signal, fftn_1);%MPO_QPC_fft = Ts.*fft(Y_est_sub(zc_tspan:end), fftn);
% signal_FFT_Mag = abs(signal_fft);
% signal_freq_compnt = [ signal_FFT_Mag(freq_pos(intrmd_dff,Fs,fftn_1),:) , signal_FFT_Mag(freq_pos(intrmd_sm,Fs,fftn_1),:) , ...
%                        signal_FFT_Mag(freq_pos(pos_freq_comp(2),Fs,fftn_1),:) ];
% spect_comp = cplng_param(1:3) ./ signal_freq_compnt;

end
%-------------------------------

%-------------------------------
function [cplng_param,model] = grd_srch_lin(signal, env, Fs, Ts, pos_freq_comp, phi, vrb, RCT)
freq_pos = @(freq,Fs,fftn) floor(freq*fftn/Fs)+1;
% freq_pos = @(freq,Fs,fftn) round(freq*fftn/Fs);
if pos_freq_comp(2) > pos_freq_comp(1)
    nb2 = [ceil( (1/pos_freq_comp(1))/(4*Ts) ) , ceil( (1/pos_freq_comp(2))/(Ts) )];
    [mods_nofrfs, model] = pac_miso_base_lin(nb2, signal, env, Fs, Ts, pos_freq_comp, phi, vrb, RCT);
    fft_mags = mods_nofrfs{1,9};
    fftn = length(fft_mags);

    intrmd_sm = pos_freq_comp(2)+pos_freq_comp(1); intrmd_dff = pos_freq_comp(2)-pos_freq_comp(1);

    low_freq_mag  = fft_mags(freq_pos(pos_freq_comp(1),Fs,fftn),:);
    high_freq_mag = fft_mags(freq_pos(pos_freq_comp(2),Fs,fftn),:);
    intrmd_sm_freq_mag  = fft_mags(freq_pos(intrmd_sm,Fs,fftn),:);
    intrmd_dff_freq_mag  = fft_mags(freq_pos(intrmd_dff,Fs,fftn),:);

	intrmd_freqs = [intrmd_sm_freq_mag,intrmd_dff_freq_mag];
	intrmd_high_freq_ratio = intrmd_freqs ./ high_freq_mag;
	modl = mean( intrmd_high_freq_ratio );
	cplng_param = [ intrmd_dff_freq_mag , intrmd_sm_freq_mag , high_freq_mag , low_freq_mag , modl ];
else
    cplng_param = [0,0,0,0,0];
    model = 0;
end
end
%-------------------------------
%%

%-------------------------------
function match_ind = intr_modul_comp(model)
mod_term_char = model{1,16}.Properties.RowNames; % Model term character strings 
n_terms = length(mod_term_char); % No. of model terms
pattern = '^u(\d+)\(t-(\d+)\)u(\d+)\(t-(\d+)\)$'; % Define the regular expression pattern
match_ind = zeros(n_terms,1);
for i = 1:n_terms
    term_str = mod_term_char(i); term_str = term_str{1}; % Term identification string
    matches = regexp(term_str, pattern, 'tokens'); % Use regexp to find matches
    % Check if matches are found and omit ux(t-c)ux(t-d)
    if ~isempty(matches)
        % Extract the integer values from the matched tokens
        I_1 = str2double(matches{1}{1});
        I_3 = str2double(matches{1}{3});
        % Check the condition I_1 ≠ I_3
        if I_1 ~= I_3
            match_ind(i) = 1;
        end
    end
end
if model{11} ~= 0 % Remove bias term index if present
    match_ind = match_ind(1:end-1);
end
end
%-------------------------------

function match_ind = lin_comp_u2(model)
mod_term_char = model{1,16}.Properties.RowNames; % Model term character strings 
n_terms = length(mod_term_char); % No. of model terms
pattern = '^u2\(t-(\d+)\)$'; % Define the regular expression pattern
match_ind = zeros(n_terms,1);
for i = 1:n_terms
    term_str = mod_term_char(i); term_str = term_str{1}; % Term identification string
    matches = regexp(term_str, pattern, 'match'); % Use regexp to find matches
    % Check if matches are found and omit ux(t-c)ux(t-d)
    if ~isempty(matches)
        match_ind(i) = 1;
    end
end
if model{11} ~= 0 % Remove bias term index if present
    match_ind = match_ind(1:end-1);
end
end

function match_ind = lin_comp_u1(model)
mod_term_char = model{1,16}.Properties.RowNames; % Model term character strings 
n_terms = length(mod_term_char); % No. of model terms
pattern = '^u1\(t-(\d+)\)$'; % Define the regular expression pattern
match_ind = zeros(n_terms,1);
for i = 1:n_terms
    term_str = mod_term_char(i); term_str = term_str{1}; % Term identification string
    matches = regexp(term_str, pattern, 'match'); % Use regexp to find matches
    % Check if matches are found and omit ux(t-c)ux(t-d)
    if ~isempty(matches)
        match_ind(i) = 1;
    end
end
if model{11} ~= 0 % Remove bias term index if present
    match_ind = match_ind(1:end-1);
end
end
%%





