function [filt_sig] = lowpass_fir(signal, Fc, Fs)
%LOWPASS_FIR Apply a zero-phase Hamming-window FIR low-pass filter.
%   A fixed 5001st-order linear-phase FIR filter is designed and applied in
%   forward and reverse directions with FILTFILT to remove phase delay.
%   Inputs
%   ------
%   signal    : Signal vector to filter.
%   Fc        : Low-pass cutoff frequency in Hz.
%   Fs        : Sampling frequency in Hz.
%   Output
%   ------
%   filt_sig  : Zero-phase low-pass-filtered signal, with the same orientation
%               as SIGNAL.

%% Design and apply the FIR filter

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
