function [s_LF, s_HF, pink, pink_aug] = pink_noise_LF_HF(N, Fs, LF_freq, HF_freq)
%PINK_NOISE_LF_HF Generate pink noise and isolated slow/fast components.
%   White noise is spectrally shaped with a 1/sqrt(f) amplitude law, then
%   zero-phase FIR bandpass filters extract the requested low- and
%   high-frequency components. This supports the non-stationary synthetic
%   oscillation construction illustrated in Figure 12 of the paper.
%   Inputs
%   ------
%   N         : Number of samples to generate.
%   Fs        : Sampling frequency in Hz.
%   LF_freq   : [lowEdge, highEdge] of the slow band in Hz.
%   HF_freq   : [lowEdge, highEdge] of the fast band in Hz.
%   Outputs
%   -------
%   s_LF      : Pink-noise component in LF_freq.
%   s_HF      : Pink-noise component in HF_freq.
%   pink      : Full spectrally shaped pink-noise realisation.
%   pink_aug  : Residual pink noise after subtracting s_LF and s_HF.

%% Generate a 1/f power-spectrum noise realisation
%----------------- Generate pink noise -----------------------
white = randn(1, N);
f = fft(white);
frequencies = [ 0:Fs/N:(Fs/2)-(Fs/N) , fliplr(0:Fs/N:(Fs/2)-(Fs/N)) ];
scaling = 1 ./ sqrt(frequencies); % 1/f amplitude decay
scaling(1)=0; scaling(end)=0;

len_f = length(f); len_scaling = length(scaling);
min_len = min([len_f, len_scaling]);
f = f(1:min_len) .* scaling(1:min_len);
pink = real(ifft(f));

%-------Low freq FIR-------
%% Extract the non-stationary low-frequency oscillation
numTaps = 5001;%501;  % Must be odd for linear phase
b = fir1(numTaps - 1, [LF_freq(1), LF_freq(2)]/(Fs/2), 'bandpass', hamming(numTaps));
s_LF = filtfilt(b, 1, pink);

%-------High freq FIR-------
%% Extract the non-stationary high-frequency oscillation
numTaps = 2501;%251;  % Shorter needed for higher freqs
b = fir1(numTaps - 1, [HF_freq(1), HF_freq(2)]/(Fs/2), 'bandpass', kaiser(numTaps, 8));
s_HF = filtfilt(b, 1, pink);


%% Form a residual noise process without the two extracted bands
pink_aug = pink - (s_LF + s_HF);
end
