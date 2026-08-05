function [data] = nrrw_bnd_fft_filt(data_raw,Fs,Fc,K,filt_typ)
N = length(data_raw);
data_fft = fft(data_raw);
w = 0:Fs/N:Fs-(Fs/N); 

% Fc_freq_ind = floor(((Fc/Fs)*N) + 1);
% Fc_freq_ind_mrr = floor((((Fs-Fc)/Fs)*N) + 1);
%freq_rmv_ind = zeros(size(w));
%freq_rmv_ind([Fc_freq_ind-K:Fc_freq_ind+K, Fc_freq_ind_mrr-K:Fc_freq_ind_mrr+K]) = 1;

Fc_freq_ind = floor( ( (Fc/Fs)*N ) )+1;
% [~, Fc_freq_ind] = min(abs(w-Fc));
f_low = w(Fc_freq_ind-K);
f_high = w(Fc_freq_ind+K);

% K = floor( ( ( (K)/Fs )*N ) + 1) * (Fs/N);
% Fc_mK_freq_ind = floor( ( ( (Fc-K)/Fs )*N ) + 1);
% Fc_pK_freq_ind = floor( ( ( (Fc+K)/Fs )*N ) + 1);
% f_low = w(Fc_mK_freq_ind);
% f_high = w(Fc_pK_freq_ind);

switch filt_typ
    case 'bw'
        W_shifted = (w>=f_low & w<=f_high) | (w>=(Fs-f_high) & w<=(Fs-f_low));
    case 'sbp'
        W_shifted = smooth_bandpass_fft_window(f_low, f_high, 0.1, N, Fs);
    case 'guss'
        W_shifted = gaussian_fft_window(f_low, f_high, w(Fc_freq_ind), N, Fs);
end

W_shifted = W_shifted./max(abs(W_shifted)); %Normalise the window such that max in 1

% figure;stem(w,W_shifted); title(num2str(Fc));

data_fft_rmv = data_fft.*W_shifted;
data = ifft(data_fft_rmv,'symmetric');
end
%%

function [W_shifted] = smooth_bandpass_fft_window(f_low, f_high, transition, fftn, Fs)

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

% function [W_shifted] = smooth_bandpass_fft_window2(f_low, f_high, Fc_freq_ind, fftn, Fs)
% 
% fft_res = Fs/fftn;
% freqs = (-Fs/2):fft_res:(Fs/2)-fft_res;
% Fc = Fc_freq_ind*fft_res;
% 
% % Define transition edges
% f1 = f_low;
% f4 = f_high;
% 
% % Initialize window
% W = zeros(size(freqs));
% 
% f = abs(freqs);
% 
% ind_transition_1 = f1 <= f & f < Fc;
% ind_passband = f==Fc;
% ind_transition_2 = Fc < f & f <= f4;
% 
% W(ind_transition_1) = 0.5 * (1 - cos(pi * (f(ind_transition_1) - f1) / (Fc - f1)));
% W(ind_passband) = 1;
% W(ind_transition_2) = 0.5 * (1 + cos(pi * (f(ind_transition_2) - Fc) / (f4 - Fc)));
% 
% % Shift for FFT compatibility
% W_shifted = ifftshift(W);  % Shift zero freq to index 1 for FFT filtering
% end


function [W_shifted] = gaussian_fft_window(f_low, f_high, f_center, fftn, Fs)

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