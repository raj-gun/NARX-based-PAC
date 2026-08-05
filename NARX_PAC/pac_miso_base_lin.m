function [mods_nofrfs, model] = pac_miso_base_lin(nb2, y_sysid, env, Fs, Ts, pos_freq_comp, phi, vrb, RCT)
%PAC_MISO_BASE_LIN Identify the linear ARX model used for grid pre-screening.
%   This is the computationally inexpensive initial-scan counterpart to
%   PAC_MISO_BASE. A two-input first-order polynomial model (nl_ord_max = 1)
%   is fitted to determine whether both the low- and high-frequency input
%   clusters are represented. The procedure implements the ARX screening step
%   in Algorithm 1 of the paper before second-order NARX identification.
%   Inputs
%   ------
%   nb2            : [maxLag_u1, maxLag_u2] in samples.
%   y_sysid        : Observed signal used as the identification output.
%   env            : N-by-2 matrix of filtered low- and high-frequency inputs.
%   Fs, Ts         : Sampling frequency in Hz and sampling interval in seconds.
%   pos_freq_comp  : [fL, fH] centre frequencies in Hz.
%   phi            : Input phases for stationary sinusoidal simulation.
%   vrb            : Verbose/plotting flag.
%   RCT            : NonSysID-i control argument passed unchanged.
%   Outputs
%   -------
%   mods_nofrfs    : Cell array of model information and simulated linear
%                    spectrum. Entries {9}, {10}, and {11} are magnitude,
%                    phase in degrees, and frequency vector, respectively.
%   model          : Identified NonSysID-i linear model or diagnostic text.
%   External dependency: NonSysID-i and its helper functions.

%% Configure stationary inputs for the ARX screening simulation
fftn_dsrd = (Fs/0.1);
N = length(y_sysid);
tspan_sysid = 0:Ts:(N*Ts-Ts);
tspan_nofrf = 0:Ts:(fftn_dsrd+100)*Ts;%-2.5:Ts:200.5;%-6+Ts:Ts:6-Ts;%tspan_sysid;%
%% Inputs

%----- filter and pass both low and high freq ----------
u = env;
u_fft = [ cos( 2.*pi.*pos_freq_comp(1).*tspan_nofrf + phi(1) )' , cos( 2.*pi.*pos_freq_comp(2).*tspan_nofrf + phi(2) )' ];
%% Normalise inputs and output

%% Scale stationary inputs to match the filtered data
std_cos = sqrt(0.5);

u_std = std(u,1); scl_fctr = std_cos./u_std;
u_fft = u_fft ./ scl_fctr;

fftn = fftn_dsrd;
w = 0:Fs/fftn:Fs-(Fs/fftn);
%% Configure and run linear NonSysID-i identification
mod_type = 'ARX-i'; % Model type ARX-i
n_inpts=2; % Specify number of inputs
na1=1;na2=1; % Maximum and minimum output lags
nb1=[1,1]; % Minimum input lags
nl_ord_max=1; % Maximum order of polynomial nonlinearity considered
x_iOFR = [false,false]; % Run more than one iteration of iOFR for [linear model ,nonlinear model]
stp_cri = {'PRESS_thresh', []}; D1_thresh = [10^(-3),[]]; % Stoping criteria for [linear model ,nonlinear model]. PRESS_thresh/BIC_thresh
is_bias=1; % Specify if bias/DC off set is required, 0, or not, 1.
KSA_h=20; % Specify the number of steps for k-steps ahead prediction 
sim=[1,1].*vrb; % Specify whether to simulate model and display results respectively
displ=0; % Set to 1 to display all models generated from iOFRs, 0 otherwise 
parall = [0,0]; % Set 1 or 0 to use parallel processing to accelerate iOFRs, for [linear model ,nonlinear model]

% Run NonSysID_i
% Failed or numerically extreme models are returned as zero-valued spectra,
% allowing the outer grid scan to treat the pair as unsuitable.
try
    [model, Mod_Val_dat, iOFR_table_lin, iOFR_table_nl, best_mod_ind_lin, best_mod_ind_nl, ~] = ...
        NonSysID_i(mod_type,u,y_sysid,na1,na2,nb1,nb2,nl_ord_max,is_bias,n_inpts,KSA_h,RCT,x_iOFR,stp_cri,D1_thresh,displ,sim,parall);
catch
    MPO_u1u2_FFT_Mag = zeros(fftn,1);
    MPO_u1u2_FFT_Phs = zeros(fftn,1);

    mods_nofrfs{1,1} = 0;
    mods_nofrfs{1,2} = 0;
    mods_nofrfs{1,3} = 0;
    mods_nofrfs{1,4} = 0;
    mods_nofrfs{1,5} = 0;
    mods_nofrfs{1,6} = 0;%G_LS_2;
    mods_nofrfs{1,7} = 0;%Y_NOFRF_MLS_n;
    mods_nofrfs{1,8} = 0;%rmv_U_mat;
    mods_nofrfs{1,9} = MPO_u1u2_FFT_Mag;
    mods_nofrfs{1,10} = MPO_u1u2_FFT_Phs;
    mods_nofrfs{1,11} = 0;

    model = 'A model cannot be built';
    return
end

if sum( abs(model{1,5}) > 1e3 ) ~= 0
    MPO_u1u2_FFT_Mag = zeros(fftn,1);
    MPO_u1u2_FFT_Phs = zeros(fftn,1);

    mods_nofrfs{1,1} = 0;
    mods_nofrfs{1,2} = 0;
    mods_nofrfs{1,3} = 0;
    mods_nofrfs{1,4} = 0;
    mods_nofrfs{1,5} = 0;
    mods_nofrfs{1,6} = 0;%G_LS_2;
    mods_nofrfs{1,7} = 0;%Y_NOFRF_MLS_n;
    mods_nofrfs{1,8} = 0;%rmv_U_mat;
    mods_nofrfs{1,9} = MPO_u1u2_FFT_Mag;
    mods_nofrfs{1,10} = MPO_u1u2_FFT_Phs;
    mods_nofrfs{1,11} = 0;

    model = 'Model parameters are too large, model is tending to be unstable';
    return
end

%% Extract the components of the model that generates quadratic phase couplings

%% Verify that both linear input clusters are present
match_ind_arx_u1 = lin_comp_u1(iOFR_table_lin,best_mod_ind_lin,model);
match_ind_arx_u2 = lin_comp_u2(iOFR_table_lin,best_mod_ind_lin,model);
match_ind_u1u2_arx = match_ind_arx_u1 | match_ind_arx_u2;

if sum(match_ind_arx_u1)~=0 && sum(match_ind_arx_u2)~=0
    [~,Y_est_sub] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_u1u2_arx);
    fftn = fftn_dsrd;
       
    MPO_u1u2_fft = Ts.*fft(Y_est_sub, fftn); %MPO_u1u2_fft = Ts.*fft(Y_est_sub(zc_tspan:end), fftn);
    MPO_u1u2_FFT_Mag = abs(MPO_u1u2_fft)./length(Y_est_sub);
    MPO_u1u2_FFT_Phs = angle(MPO_u1u2_fft).*(180/pi);
else
    MPO_u1u2_FFT_Mag = zeros(fftn,1);
    MPO_u1u2_FFT_Phs = zeros(fftn,1);
end

%% Optional diagnostic plots and model display
if vrb==1
    figure;ax=subplot(2,1,1);plot(w, MPO_u1u2_FFT_Mag,  'k', 'LineWidth', 1.5);axis([0,Fs/2,0,inf]);
    set(ax,'YTick',[]); box(ax,'off');
    subplot(2,1,2);plot(w, MPO_u1u2_FFT_Phs);axis([0,Fs/2,0,inf]);
    disp(iOFR_table_lin{best_mod_ind_lin,1});
end
%% Package validation statistics and simulated spectrum
mods_nofrfs = cell(1,8);
if size(iOFR_table_nl)~=[1,1]
    mod_info = iOFR_table_nl{best_mod_ind_nl,1};
    mod_val_temp = iOFR_table_nl{best_mod_ind_nl,16}([2,4,5],:);
    mod_val_nl = sum(mod_val_temp(2:end,:),1);
    mod_val_lin = mod_val_temp(1,:);
    msse = model{14};
    model{1,length(model)+1} = iOFR_table_nl{best_mod_ind_nl,1}.Properties.RowNames;
else
    val_stats = mod_val_stats(Mod_Val_dat);
    mod_info = 0; mod_val_nl = [1e10 1e10 1e10]; mod_val_lin = val_stats(2,:); msse = model{14};
end
mods_nofrfs{1,1} = mod_info;
mods_nofrfs{1,2} = msse;
mods_nofrfs{1,3} = mod_val_nl;
mods_nofrfs{1,4} = mod_val_lin;
mods_nofrfs{1,5} = sse;
mods_nofrfs{1,6} = 0;%G_LS_2;
mods_nofrfs{1,7} = 0;%Y_NOFRF_MLS_n;
mods_nofrfs{1,8} = 0;%rmv_U_mat;
mods_nofrfs{1,9} = MPO_u1u2_FFT_Mag;
mods_nofrfs{1,10} = MPO_u1u2_FFT_Phs;
mods_nofrfs{1,11} = w;
end

%% Local functions
function match_ind = intr_modul_comp(iOFR_table_nl,best_mod_ind_nl,model)
%INTR_MODUL_COMP Select cross-input quadratic terms from a model table.
%   Inputs identify the selected NonSysID-i model table and bias convention.
%   Output MATCH_IND selects cross-input quadratic terms.
%   Retained for compatibility with related NARX routines; the linear scan
%   does not require this selector.
mod_term_char = iOFR_table_nl{best_mod_ind_nl,1}.Properties.RowNames; % Model term character strings 
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

function match_ind = lin_comp(iOFR_table_nl,best_mod_ind_nl,model)
%LIN_COMP Select all first-order delayed-input terms from a model table.
%   Inputs identify the selected model table. Output MATCH_IND selects all
%   first-order delayed-input terms.
mod_term_char = iOFR_table_nl{best_mod_ind_nl,1}.Properties.RowNames; % Model term character strings 
n_terms = length(mod_term_char); % No. of model terms
pattern = '^u(\d+)\(t-(\d+)\)$'; % Define the regular expression pattern
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

function match_ind = lin_comp_u2(iOFR_table_nl,best_mod_ind_nl,model)
%LIN_COMP_U2 Select high-frequency input terms Sigma_u2.
%   Inputs identify the selected model table. Output MATCH_IND selects the
%   high-frequency linear-input cluster.
mod_term_char = iOFR_table_nl{best_mod_ind_nl,1}.Properties.RowNames; % Model term character strings 
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

function match_ind = lin_comp_u1(iOFR_table_nl,best_mod_ind_nl,model)
%LIN_COMP_U1 Select low-frequency input terms Sigma_u1.
%   Inputs identify the selected model table. Output MATCH_IND selects the
%   low-frequency linear-input cluster.
mod_term_char = iOFR_table_nl{best_mod_ind_nl,1}.Properties.RowNames; % Model term character strings 
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
