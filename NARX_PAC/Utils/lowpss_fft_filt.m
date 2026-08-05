function [data] = lowpss_fft_filt(data_raw,Fs,Fct,transition_gap)
N = length(data_raw);
data_fft = fft(data_raw);
w = 0:Fs/N:Fs-(Fs/N); %freq_rmv_ind = (w>=f1 & w<=f2) | (w>=(Fs-f2) & w<=(Fs-f1));

% Fc_freq_ind = floor(((Fc/Fs)*N) + 1);
% Fc_freq_ind_mrr = floor((((Fs-Fc)/Fs)*N) + 1);
%freq_rmv_ind = zeros(size(w));
%freq_rmv_ind([Fc_freq_ind-K:Fc_freq_ind+K, Fc_freq_ind_mrr-K:Fc_freq_ind_mrr+K]) = 1;

% [~, Fc_freq_ind] = min(abs(w-Fct));
Fc_freq_ind = floor( ( (Fct/Fs)*N ) )+1;

W_shifted = smooth_lwpss_fft_window(w(Fc_freq_ind), transition_gap, N, Fs); 

W_shifted = W_shifted./max(abs(W_shifted)); %Normalise the window such that max in 1

% figure;stem(w,W_shifted);

data_fft_rmv = data_fft.*W_shifted;
data = ifft(data_fft_rmv,'symmetric');
end

function [W_shifted] = smooth_lwpss_fft_window(Fct, transition_gap, fftn, Fs)

fft_res = Fs/fftn;
freqs = (-Fs/2):fft_res:(Fs/2)-fft_res;

% Transition width
t_width = transition_gap;
t_width = round(t_width/fft_res) * fft_res; if t_width==0; t_width=fft_res; end

% Define transition edges
f3 = Fct - t_width;
f4 = Fct + t_width;

% Initialize window
W = zeros(size(freqs));

f = abs(freqs);

ind_passband = f <= Fct;
ind_transition_2 = f3 < f & f <= f4;

W(ind_passband) = 1;
W(ind_transition_2) = 0.5 * (1 + cos(pi * (f(ind_transition_2) - f3) / (f4 - f3)));

% Shift for FFT compatibility
W_shifted = ifftshift(W);  % Shift zero freq to index 1 for FFT filtering
end