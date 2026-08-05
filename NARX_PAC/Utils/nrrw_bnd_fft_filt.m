function [data] = nrrw_bnd_fft_filt(data_raw,Fs,Fc,K,filt_typ)
%NRRW_BND_FFT_FILT Extract a narrow frequency band in the FFT domain.
%   A real, two-sided spectral window is centred at Fc and applied directly
%   to the FFT. The even window and symmetric inverse FFT preserve zero phase,
%   supporting the practical input extraction in Section III F of the paper.
%   Inputs
%   ------
%   data_raw  : Signal vector to filter.
%   Fs        : Sampling frequency in Hz.
%   Fc        : Centre frequency in Hz.
%   K         : Half-bandwidth in FFT bins; passband edges are Fc +/- K bins.
%   filt_typ  : 'bw' for a brick-wall window, 'sbp' for a 10% cosine-roll-off
%               smooth bandpass, or 'guss' for a Gaussian window.
%   Output
%   ------
%   data      : Real narrowband signal obtained by symmetric inverse FFT.

%% Determine the discrete passband edges and construct the selected window
N = length(data_raw);
data_fft = fft(data_raw);
w = 0:Fs/N:Fs-(Fs/N); 


Fc_freq_ind = floor( ( (Fc/Fs)*N ) )+1;
f_low = w(Fc_freq_ind-K);
f_high = w(Fc_freq_ind+K);


switch filt_typ
    case 'bw'
        W_shifted = (w>=f_low & w<=f_high) | (w>=(Fs-f_high) & w<=(Fs-f_low));
    case 'sbp'
        W_shifted = smooth_bandpass_fft_window(f_low, f_high, 0.1, N, Fs);
    case 'guss'
        W_shifted = gaussian_fft_window(f_low, f_high, w(Fc_freq_ind), N, Fs);
end

W_shifted = W_shifted./max(abs(W_shifted)); %Normalise the window such that max in 1


%% Apply the two-sided window and return the real filtered signal
data_fft_rmv = data_fft.*W_shifted;
data = ifft(data_fft_rmv,'symmetric');
end

function [W_shifted] = smooth_bandpass_fft_window(f_low, f_high, transition, fftn, Fs)
%SMOOTH_BANDPASS_FFT_WINDOW Build an even cosine-roll-off bandpass window.
%   Inputs define the passband edges, fractional transition width, FFT
%   length, and sampling frequency. Output W_SHIFTED follows MATLAB FFT order.
%   TRANSITION is expressed as a fraction of the requested passband width.

fft_res = Fs/fftn;
freqs = (-Fs/2):fft_res:(Fs/2)-fft_res;

% Transition width
bw = f_high - f_low;
t_width = transition * bw;
t_width = round(t_width/fft_res) * fft_res; if t_width==0; t_width=fft_res; end

% Define transition edges
f1 = f_low - t_width;
f2 = f_low + t_width;
f3 = f_high - t_width;
f4 = f_high + t_width;

% Initialize window
W = zeros(size(freqs));

f = abs(freqs);

ind_transition_1 = f1 <= f & f < f2;
ind_passband = f2 <= f & f <= f3;
ind_transition_2 = f3 < f & f <= f4;

W(ind_transition_1) = 0.5 * (1 - cos(pi * (f(ind_transition_1) - f1) / (f2 - f1)));
W(ind_passband) = 1;
W(ind_transition_2) = 0.5 * (1 + cos(pi * (f(ind_transition_2) - f3) / (f4 - f3)));

% Shift for FFT compatibility
W_shifted = ifftshift(W);  % Shift zero freq to index 1 for FFT filtering
end


function [W_shifted] = gaussian_fft_window(f_low, f_high, f_center, fftn, Fs)
%GAUSSIAN_FFT_WINDOW Build an even Gaussian bandpass window.
%   Inputs define the nominal band, centre frequency, FFT length, and
%   sampling rate. Output W_SHIFTED follows MATLAB FFT order.
%   The requested band width is interpreted as the Gaussian full width at
%   half maximum and converted to the corresponding standard deviation.

fft_res = Fs/fftn;
freqs = (-Fs/2):Fs/fftn:(Fs/2)-(Fs/fftn);

f_bw = f_high - f_low; % Frequency bandwidth
f_bw = round(f_bw/fft_res) * fft_res; if f_bw==0; f_bw=fft_res; end % Adjust bandwidth to nearest frequency bin gap
% Calculate Gaussian width parameter (standard deviation)
sigma = f_bw / (2*sqrt(2*log(2)));  % Convert FWHM to sigma

% Two-sided Gaussian window
W = exp(-0.5 * ((abs(freqs) - f_center)/sigma).^2);

% Shift for FFT compatibility
W_shifted = ifftshift(W);  % Shift zero freq to index 1 for FFT filtering
end
