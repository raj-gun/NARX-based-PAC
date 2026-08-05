function [data] = lowpss_fft_filt(data_raw,Fs,Fct,transition_gap)
%LOWPSS_FFT_FILT Apply a smooth zero-phase low-pass filter in the FFT domain.
%   The routine multiplies the two-sided spectrum by an even cosine-roll-off
%   window and returns the real inverse FFT. The symmetric frequency response
%   avoids an imposed phase shift.
%   Inputs
%   ------
%   data_raw        : Signal vector to filter.
%   Fs              : Sampling frequency in Hz.
%   Fct             : Passband cutoff frequency in Hz.
%   transition_gap  : Transition half-width in Hz around Fct.
%   Output
%   ------
%   data            : Low-pass-filtered real signal.

%% Construct the low-pass window and filter the spectrum
N = length(data_raw);
data_fft = fft(data_raw);
w = 0:Fs/N:Fs-(Fs/N); %freq_rmv_ind = (w>=f1 & w<=f2) | (w>=(Fs-f2) & w<=(Fs-f1));

Fc_freq_ind = floor( ( (Fct/Fs)*N ) )+1;

W_shifted = smooth_lwpss_fft_window(w(Fc_freq_ind), transition_gap, N, Fs); 

W_shifted = W_shifted./max(abs(W_shifted)); %Normalise the window such that max in 1


data_fft_rmv = data_fft.*W_shifted;
data = ifft(data_fft_rmv,'symmetric');
end

function [W_shifted] = smooth_lwpss_fft_window(Fct, transition_gap, fftn, Fs)
%SMOOTH_LWPSS_FFT_WINDOW Build an even cosine-roll-off low-pass window.
%   Inputs define the cutoff, transition width, FFT length, and sampling
%   rate. Output W_SHIFTED follows MATLAB FFT order.
%   The returned vector is shifted to MATLAB FFT ordering, with zero
%   frequency at index 1.

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
