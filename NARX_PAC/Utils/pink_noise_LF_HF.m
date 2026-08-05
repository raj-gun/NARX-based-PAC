function [s_LF, s_HF, pink, pink_aug] = pink_noise_LF_HF(N, Fs, LF_freq, HF_freq)
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
%-------------------------------------------------------------

%-------Low freq FIR-------
numTaps = 5001;%501;  % Must be odd for linear phase
b = fir1(numTaps - 1, [LF_freq(1), LF_freq(2)]/(Fs/2), 'bandpass', hamming(numTaps));
s_LF = filtfilt(b, 1, pink);
%--------------------------

%-------High freq FIR-------
numTaps = 2501;%251;  % Shorter needed for higher freqs
b = fir1(numTaps - 1, [HF_freq(1), HF_freq(2)]/(Fs/2), 'bandpass', kaiser(numTaps, 8));
s_HF = filtfilt(b, 1, pink);
%---------------------------

% % -------Low freq FFT-------
% s_LF = fft_bndpss_flt(pink, Fs, LF_freq(1), LF_freq(2));
% % -------High freq FFT-------
% s_HF = fft_bndpss_flt(pink, Fs, HF_freq(1), HF_freq(2));
% % ---------------------------

pink_aug = pink - (s_LF + s_HF);
end