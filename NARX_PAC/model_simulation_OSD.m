function [LF_sig, AM_sig] = model_simulation_OSD(model,pos_freq_comp,Ts,phi,scl_fctr)
Fs = 1/Ts;
fftn_dsrd = (Fs/0.1);
tspan_nofrf = 0:Ts:(fftn_dsrd+100)*Ts;%-2.5:Ts:200.5;%-6+Ts:Ts:6-Ts;%tspan_sysid;%
% tspan_nofrf = -2.5:Ts:200.5;%-6+Ts:Ts:6-Ts;%tspan_sysid;%
u_fft = [ scl_fctr(1).*cos( 2.*pi.*pos_freq_comp(1).*tspan_nofrf + phi(1) )' , scl_fctr(2).*cos( 2.*pi.*pos_freq_comp(2).*tspan_nofrf + phi(2) )' ];

match_ind_intmod = intr_modul_comp(model);
match_ind_arx_u2 = lin_comp_u2(model);
match_ind_AM_narx = match_ind_intmod | match_ind_arx_u2;
match_ind_LF_arx = lin_comp_u1(model);

[~,LF_sig] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_LF_arx); % Low-freq signal
% tspan_nofrf_trim = tspan_nofrf( length(u_fft) - length(LF_sig) + 1:end );
% [~,zc_tspan] = min(abs(tspan_nofrf_trim));
% LF_sig = LF_sig(zc_tspan:end,1); 

[~,AM_sig] = model_simulation_clstr(model,u_fft,u_fft.*0,match_ind_AM_narx); % High-freq amplitude modulated signal
% AM_sig = AM_sig(zc_tspan:end,1);
end

function match_ind = intr_modul_comp(model)
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