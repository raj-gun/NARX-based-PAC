function [mods_nofrfs, model] = pac_miso_base(nb2, y_sysid, env, Fs, Ts, pos_freq_comp, phi, vrb, RCT)
%PAC_MISO_BASE Identify the second-order canonical NARX model for one pair.
%   A two-input, input-only polynomial NARX model is identified from the
%   band-limited low- and high-frequency inputs in ENV and the observed signal
%   Y_SYSID. The identified model is then simulated with stationary sinusoidal
%   inputs. Only Sigma_u1, Sigma_u2, and Sigma_u1u2 are retained to form the
%   canonical PAC approximation described in Sections III B-C and Figure 10.
%   Inputs
%   ------
%   nb2            : [maxLag_u1, maxLag_u2] in samples. The caller normally
%                    uses one quarter slow period and one full fast period.
%   y_sysid        : Observed signal used as the identification output.
%   env            : N-by-2 matrix of filtered low- and high-frequency inputs.
%   Fs, Ts         : Sampling frequency in Hz and sampling interval in seconds.
%   pos_freq_comp  : [fL, fH] centre frequencies in Hz.
%   phi            : [phiL, phiH] phases for canonical sinusoidal simulation.
%   vrb            : Verbose/plotting flag passed to validation and diagnostics.
%   RCT            : NonSysID-i control argument passed unchanged.
%   Outputs
%   -------
%   mods_nofrfs    : Cell array containing model/validation information and
%                    canonical simulation results. Key entries are:
%                    {9} magnitude spectrum, {10} phase spectrum in degrees,
%                    {11} frequency vector, and {12} simulated columns
%                    [Sigma_u1, Sigma_u2, Sigma_u1u2].
%   model          : Identified NonSysID-i model cell array, or a diagnostic
%                    text message when a valid nonlinear model is unavailable.
%   External dependency: NonSysID-i and its helper functions.

%% Configure the stationary canonical simulation
fftn_dsrd = (Fs/0.1);
N = length(y_sysid);
tspan_sysid = 0:Ts:(N*Ts-Ts);
tspan_nofrf = 0:Ts:(fftn_dsrd+1000)*Ts;%-2.5:Ts:200.5;%-6+Ts:Ts:6-Ts;%tspan_sysid;%
%% Inputs

%----- filter and pass both low and high freq ----------
u = env;
u_fft = [ cos( 2.*pi.*pos_freq_comp(1).*tspan_nofrf + phi(1) )' , cos( 2.*pi.*pos_freq_comp(2).*tspan_nofrf + phi(2) )' ];
%% Normalise inputs and output

%% Match canonical input amplitudes to the filtered identification inputs
% A unit cosine has standard deviation sqrt(1/2); rescaling preserves the
% relative amplitudes specified in Section III F of the paper.
std_cos = sqrt(0.5);


u_std = std(u,1); scl_fctr = std_cos./u_std;
u_fft = u_fft ./ scl_fctr;


fftn = fftn_dsrd;
w = 0:Fs/fftn:Fs-(Fs/fftn);
%% Configure and run second-order NonSysID-i identification
mod_type = 'ARX-i'; % Model type ARX-i
n_inpts=2; % Specify number of inputs
na1=1;na2=1; % Maximum and minimum output lags
nb1=[1,1]; % Minimum input lags
nl_ord_max=2; % Maximum order of polynomial nonlinearity considered
x_iOFR = [false,false]; % Run more than one iteration of iOFR for [linear model ,nonlinear model]
stp_cri = {'PRESS_thresh', 'PRESS_min'}; D1_thresh = [10^(-3),[]]; %[10^(-3),10^(-6)]; % Stoping criteria for [linear model ,nonlinear model]. PRESS_thresh/BIC_thresh
is_bias=1; % Specify if bias/DC off set is required, 0, or not, 1.
KSA_h=20; % Specify the number of steps for k-steps ahead prediction
sim=[1,1].*vrb; % Specify whether to simulate model and display results respectively
displ=0; % Set to 1 to display all models generated from iOFRs, 0 otherwise
parall = [0,0]; % Set 1 or 0 to use parallel processing to accelerate iOFRs, for [linear model ,nonlinear model]

% Run NonSysID_i
% Identification failures are converted to zero-valued spectral outputs so
% the surrounding frequency-grid search can continue safely.
try
    [model, Mod_Val_dat, iOFR_table_lin, iOFR_table_nl, best_mod_ind_lin, best_mod_ind_nl, ~] = ...
        NonSysID_i(mod_type,u,y_sysid,na1,na2,nb1,nb2,nl_ord_max,is_bias,n_inpts,KSA_h,RCT,x_iOFR,stp_cri,D1_thresh,displ,sim,parall);
catch
    MPO_QPC_FFT_Mag = zeros(fftn,1);
    MPO_QPC_FFT_Phs = zeros(fftn,1);

    mods_nofrfs{1,1} = 0;
    mods_nofrfs{1,2} = 0;
    mods_nofrfs{1,3} = 0;
    mods_nofrfs{1,4} = 0;
    mods_nofrfs{1,5} = 0;
    mods_nofrfs{1,6} = 0;%G_LS_2;
    mods_nofrfs{1,7} = 0;%Y_NOFRF_MLS_n;
    mods_nofrfs{1,8} = 0;%rmv_U_mat;
    mods_nofrfs{1,9} = MPO_QPC_FFT_Mag;
    mods_nofrfs{1,10} = MPO_QPC_FFT_Phs;
    mods_nofrfs{1,11} = 0;
    mods_nofrfs{1,12} = 0;

    model = 'A model cannot be built';
    return
end

%% Reject models that do not provide a usable nonlinear canonical structure
if size(iOFR_table_nl)==[1,1]
    MPO_QPC_FFT_Mag = zeros(fftn,1);
    MPO_QPC_FFT_Phs = zeros(fftn,1);

    mods_nofrfs{1,1} = 0;
    mods_nofrfs{1,2} = 0;
    mods_nofrfs{1,3} = 0;
    mods_nofrfs{1,4} = 0;
    mods_nofrfs{1,5} = 0;
    mods_nofrfs{1,6} = 0;%G_LS_2;
    mods_nofrfs{1,7} = 0;%Y_NOFRF_MLS_n;
    mods_nofrfs{1,8} = 0;%rmv_U_mat;
    mods_nofrfs{1,9} = MPO_QPC_FFT_Mag;
    mods_nofrfs{1,10} = MPO_QPC_FFT_Phs;
    mods_nofrfs{1,11} = 0;
    mods_nofrfs{1,12} = 0;

    model = 'Model is linear';
    return
end

if sum( abs(model{1,5}) > 1e3 ) ~= 0
    MPO_QPC_FFT_Mag = zeros(fftn,1);
    MPO_QPC_FFT_Phs = zeros(fftn,1);

    mods_nofrfs{1,1} = 0;
    mods_nofrfs{1,2} = 0;
    mods_nofrfs{1,3} = 0;
    mods_nofrfs{1,4} = 0;
    mods_nofrfs{1,5} = 0;
    mods_nofrfs{1,6} = 0;%G_LS_2;
    mods_nofrfs{1,7} = 0;%Y_NOFRF_MLS_n;
    mods_nofrfs{1,8} = 0;%rmv_U_mat;
    mods_nofrfs{1,9} = MPO_QPC_FFT_Mag;
    mods_nofrfs{1,10} = MPO_QPC_FFT_Phs;
    mods_nofrfs{1,11} = 0;
    mods_nofrfs{1,12} = 0;

    model = 'Model parameters are too large, model is tending to be unstable';
    return
end

%% Extract the components of the model that generates quadratic phase couplings

%% Locate Sigma_u1, Sigma_u2, and Sigma_u1u2 term clusters
match_ind_u1 = lin_comp_u1(model);
match_ind_u2 = lin_comp_u2(model);

match_ind_intmod = intr_modul_comp(model);
match_ind_arx_u1 = lin_comp_u1(model); 
match_ind_arx_u2 = lin_comp_u2(model);
match_ind_arx = lin_comp(model);
match_ind_QPC_narx = match_ind_intmod | match_ind_arx;

%% Simulate the canonical cluster responses and compute their spectrum
if sum(match_ind_intmod)~=0 && sum(match_ind_arx_u1)~=0 && sum(match_ind_arx_u2)~=0
    [~,Y_est_sub] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_QPC_narx);
    [~,Y_u1] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_u1);
    [~,Y_u2] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_u2);
    [~,Y_intmd] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_intmod);
    [~,loc] = findpeaks(Y_u1); loc = loc(1);
    Y_est_sub = Y_est_sub(loc:end , :);
    Y_u1      = Y_u1(loc:end , :);
    Y_u2      = Y_u2(loc:end , :);
    Y_intmd   = Y_intmd(loc:end , :);
    Y_u1_u2_phs_dat = [ Y_u1 , Y_u2 , Y_intmd];
    Y_u1_u2_phs_dat = Y_u1_u2_phs_dat(1:fftn_dsrd,:);
    fftn = fftn_dsrd;
    MPO_QPC_fft = (Ts/length(Y_est_sub)).*fft(Y_est_sub, fftn);%MPO_QPC_fft = Ts.*fft(Y_est_sub(zc_tspan:end), fftn);
    MPO_QPC_FFT_Mag = abs(MPO_QPC_fft);
    MPO_QPC_FFT_Phs = angle(MPO_QPC_fft).*(180/pi);
else
    MPO_QPC_FFT_Mag = zeros(fftn,1);
    MPO_QPC_FFT_Phs = zeros(fftn,1);
    Y_u1_u2_phs_dat = zeros(N,5);
end

%% Optional diagnostic plots and model tables
if vrb==1
    figure;ax=subplot(2,1,1);plot(w, MPO_QPC_FFT_Mag,  'k', 'LineWidth', 1.5);axis([0,Fs/2,0,inf]);
    set(ax,'YTick',[]); box(ax,'off');
    subplot(2,1,2);plot(w, MPO_QPC_FFT_Phs);axis([0,Fs/2,0,inf]);
    disp(iOFR_table_lin{best_mod_ind_lin,1});
    if size(iOFR_table_nl)~=[1,1]
        disp(iOFR_table_nl{best_mod_ind_nl,1});
        disp(iOFR_table_nl{best_mod_ind_nl,16}([2,4,5],:));
        disp(sum(iOFR_table_nl{best_mod_ind_nl,16}([2,4,5],:)));
    end
end
%% Package validation statistics and canonical simulation outputs
mods_nofrfs = cell(1,12);
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
mods_nofrfs{1,5} = 0;
mods_nofrfs{1,6} = 0;%G_LS_2;
mods_nofrfs{1,7} = 0;%Y_NOFRF_MLS_n;
mods_nofrfs{1,8} = 0;%rmv_U_mat;
mods_nofrfs{1,9} = MPO_QPC_FFT_Mag;
mods_nofrfs{1,10} = MPO_QPC_FFT_Phs;
mods_nofrfs{1,11} = w;
mods_nofrfs{1,12} = Y_u1_u2_phs_dat;
end

%% Local functions
function match_ind = intr_modul_comp(model)
%INTR_MODUL_COMP Select cross-input quadratic terms Sigma_u1u2.
%   Input: MODEL is a NonSysID-i model cell array. Output: MATCH_IND is a
%   logical-compatible selector aligned with the final coefficient vector.
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

function match_ind = lin_comp(model)
%LIN_COMP Select all first-order delayed-input terms.
%   Input: MODEL is a NonSysID-i model cell array. Output: MATCH_IND selects
%   every first-order delayed-input term.
mod_term_char = model{1,16}.Properties.RowNames; % Model term character strings
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

function match_ind = lin_comp_u2(model)
%LIN_COMP_U2 Select first-order terms from the high-frequency input u2.
%   Input: MODEL is a NonSysID-i model cell array. Output: MATCH_IND selects
%   only first-order delayed u2 terms.
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
%LIN_COMP_U1 Select first-order terms from the low-frequency input u1.
%   Input: MODEL is a NonSysID-i model cell array. Output: MATCH_IND selects
%   only first-order delayed u1 terms.
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
