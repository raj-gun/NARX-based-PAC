function [s_final, snr] = synthesize_pac(noise_lev)

% adapted from Kramer et al. (2008), Jrn. Nrsc. Methds.  

% this functions simulates a signal containing phase-amplitude cooupling (PAC)
% between the phase of a low frequency (15-16 Hz) and amplitude of a
% high-freqency band (60-90 Hz) for a chosen noise level. Hamming tapered
% high frequency signal is added at each cycle of low frequency component.
% sampling frequency is chosen as 1000 Hz.

% additionally 5 separate sinusoids (having NO coupling) to evaluate PAC
% estimation methods' discrimination capability

% noise lev : parameter describing the noise power
% snr : signal-to-noise ratio
% s_final: the synthesized signal

% For details of the simulation, please see (Ozkurt and Schnitzler, 2011)

dt = 0.001;
s = [];
for k=1:400
    f = rand()*1.0+15.0;                 %Create the low freq (15 Hz) signal.
    s1 = cos(2.0*pi*(0:dt*f:1-dt*f));
    good = find(s1 < -0.99);
    s2 = zeros(1,length(s1));
 
    stemp = randn(1,3000);                 %Create noisy data.
    stemp = eegfilt(stemp, 1000, 60, 90);%Make high freq (50-80 Hz) signal.
    stemp = stemp(2000:2039);             %Duration 50 ms.
    stemp = 5*hanning(40)'.*stemp;      %Hanning tapered.
    
    rindex = ceil(rand()*2);       %Add the high frequency in,
						%when the low frequency is near 0 phase.
    s2(rindex + good(1) - 20:rindex+good(1)+20-1)=stemp;
    s = [s, s1+1*s2];
end



n = noise_lev*randn(1,length(s));
snr = 10 * log10((s*s') / (n*n'))

s = s + n;
sorig = s;
s = s(1:12000);
sCropped = s(1000:11000-1);

% here we add sinusoids that do not have coupling
t = 1:10000;
fs = 1000;
sig1 = 2*sin(2*pi*t*20/fs); sig2 = 2*sin(2*pi*t*25/fs);
sig3 = 3*sin(2*pi*t*30/fs); sig4 = 3*sin(2*pi*t*40/fs);
sig5 = 4*sin(2*pi*t*100/fs);

s_final = sCropped + sig1 + sig2 + sig3 + sig4 + sig5;