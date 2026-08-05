function [LF_sig, AM_sig] = model_simulation_OSD(model,pos_freq_comp,Ts,phi,scl_fctr)
%MODEL_SIMULATION_OSD Decompose a canonical NARX-PAC model into two signals.
%   Stationary sinusoidal inputs are applied to the identified model. The
%   Sigma_u1 cluster is simulated as the low-frequency component, while the
%   union of Sigma_u2 and Sigma_u1u2 is simulated as the amplitude-modulated
%   high-frequency component. This is the decomposition illustrated in
%   Figure 10 and described in Section III C of the paper.
%   Inputs
%   ------
%   model          : Identified NonSysID-i NARX model.
%   pos_freq_comp  : [fL, fH] centre frequencies in Hz.
%   Ts             : Sampling interval in seconds.
%   phi            : [phiL, phiH] input phases in radians.
%   scl_fctr       : [AL, AH] amplitudes applied to the stationary sinusoidal
%                    inputs, normally derived from filtered-signal standard
%                    deviations as described in Section III F.
%   Outputs
%   -------
%   LF_sig         : Simulated Sigma_u1 low-frequency component.
%   AM_sig         : Simulated Sigma_u2 + Sigma_u1u2 high-frequency component.

%% Generate stationary low- and high-frequency inputs
Fs = 1/Ts;
fftn_dsrd = (Fs/0.1);
tspan_nofrf = 0:Ts:(fftn_dsrd+100)*Ts;%-2.5:Ts:200.5;%-6+Ts:Ts:6-Ts;%tspan_sysid;%
u_fft = [ scl_fctr(1).*cos( 2.*pi.*pos_freq_comp(1).*tspan_nofrf + phi(1) )' , scl_fctr(2).*cos( 2.*pi.*pos_freq_comp(2).*tspan_nofrf + phi(2) )' ];

%% Identify the canonical term clusters from their symbolic model terms
match_ind_intmod = intr_modul_comp(model);
match_ind_arx_u2 = lin_comp_u2(model);
match_ind_AM_narx = match_ind_intmod | match_ind_arx_u2;
match_ind_LF_arx = lin_comp_u1(model);

%% Simulate the low-frequency and amplitude-modulated components
[~,LF_sig] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_LF_arx); % Low-freq signal

[~,AM_sig] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_AM_narx); % High-freq amplitude modulated signal
end

function match_ind = intr_modul_comp(model)
%INTR_MODUL_COMP Select cross-input quadratic terms Sigma_u1u2.
%   Input: MODEL is an identified model with symbolic term names in its final
%   cell entry. Output MATCH_IND selects cross-input quadratic terms.
%   Terms must contain one delayed u1 factor and one delayed u2 factor;
%   same-input quadratic terms are excluded.
mod_term_char = model{1,end}; % Model term character strings 
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

function match_ind = lin_comp_u2(model)
%LIN_COMP_U2 Select linear delayed terms belonging to Sigma_u2.
%   Input: MODEL is an identified model. Output MATCH_IND selects linear u2
%   terms in the high-frequency component.
mod_term_char = model{1,end}; % Model term character strings 
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
%LIN_COMP_U1 Select linear delayed terms belonging to Sigma_u1.
%   Input: MODEL is an identified model. Output MATCH_IND selects linear u1
%   terms in the low-frequency component.
mod_term_char = model{1,end}; % Model term character strings 
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
