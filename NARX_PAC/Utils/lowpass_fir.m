function [filt_sig] = lowpass_fir(signal, Fc, Fs)

% Parameters
% Fs - Sampling frequency (Hz)
% Fc - Cutoff frequency (Hz)
N = 5001;          % Filter order (number of taps - 1)
                  % N should be even for symmetric FIR (linear phase)

% Normalized cutoff frequency (0 to 1, where 1 = Nyquist)
Wn = Fc / (Fs / 2);

% Design FIR filter using Hamming window
b = fir1(N, Wn, 'low', hamming(N+1));

filt_sig = filtfilt(b, 1, signal);

end