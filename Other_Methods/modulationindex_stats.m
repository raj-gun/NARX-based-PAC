%   adapted from the code written for Kramer et al. (2008) 
%    by Tolga Ozkurt (2011)

%   in order to implement statistically normalized modulation index (MI)
%   estimate (Canolty, 2006)
%  200 surrogate data are produced
%
%   INPUTS:
%    d      : The unfiltered raw data.
%    srate  : The sampling rate (e.g., 1000 Hz)
%    flow1  : chosen low limits of the bands for low-frequency component,
%    e.g., [flow1 flow1+bwlow]
%    fhigh1 : chosen low limits of the bands for high-frequency component,
%    e.g., [fhigh1 fhigh1+bwhigh]
%
%    bwlow, bwhigh: bandwidhts for the bandpass FIR filters
%
%   OUTPUS:
%    
%    mod2d_norm  = 2D statistically normalized modulation index (MI)
%    estimate (Canolty, 2006)
%
%
%    flow   = The frequency axis for the lower phase frequencies.
%    fhigh  = The frequency axis for the higher amplitude frequencies.
 
function [mod2d_norm, p_val_mtx, flow, fhigh] = modulationindex_stats(d, srate,flow1,fhigh1,bwlow,bwhigh)
 
 
  flow2 = flow1+bwlow;                                %2-41 Hz with 1 Hz steps.
  fhigh2 = fhigh1+bwhigh;                            %5-205 Hz with 5 Hz steps.
  
  mod2d = zeros(length(flow1), length(fhigh1));
  p_val_mtx = zeros(length(flow1), length(fhigh1)); % I added this part of p_values!
  
  for i=1:length(flow1)
      theta=eegfilt(d,srate,flow1(i),flow2(i));  %Compute the low freq signal.
      theta=theta(srate:length(d)-srate-1);           %Drop the first and last second.
      phase = angle(hilbert(theta));                  %Compute the low freq phase.
      ['Loops remaining = ', num2str(length(flow1)-i)]
      for j=1:length(fhigh1)
        gamma=eegfilt(d,srate,fhigh1(j),fhigh2(j));%Compute the high freq signal.
        gamma=gamma(srate:length(d)-srate-1);         %Drop the first and last second.
        amp = abs(hilbert(gamma));                 %Compute the high freq amplitude.

        %Compute the modulation index.
        [m_norm_length, p_val] = compute_modulation_index_norm(amp, phase, srate);
        mod2d(i,j) = m_norm_length;
        p_val_mtx(i,j) = p_val;
      end
  end
 
  
  
  mod2d_norm = mod2d; % here I threshold the ones having no significance
  mod2d_norm(p_val_mtx > 0.01) = 0;
  
  %Plot the thresholded two-dimensional modulation index.
  figure
  flow = (flow1 + flow2) / 2;
  fhigh = (fhigh1 + fhigh2) / 2;
  imagesc(flow, fhigh, mod2d_norm');  colorbar;
  axis xy
  set(gca, 'FontSize', 18);
  xlabel('Phase Frequency [Hz]');  ylabel('Amplitude Frequency [Hz]');
  
  %% end of the function
  
 

 
function [m_norm_length, p_values, m_raw] = compute_modulation_index_norm(a, p, srate)
 
  numpoints = length(a);
  numsurrogate=200;                                 %Use 200 surrogates.
  
  minskip=srate;                                  %The min amount to shift amplitude.
  maxskip=numpoints-srate;                        %The max amount to shift amplitude.
  skip = ceil(numpoints.*rand(numsurrogate*2,1)); %Create a list of amplitude shifts.
  skip(find(skip > maxskip))=[];                  %Make sure they're not too big.
  skip(find(skip < minskip))=[];                  %Or too small.
  skip = skip(1:numsurrogate,1);
  
  z = a.*exp(1i*p);
  m_raw = mean(z);   %Compute the mean length of unshuffled data.
  

 
  surrogate_m = zeros(numsurrogate, 1);
  for s=1:numsurrogate                                      %For each surrogate,
   surrogate_amplitude = [a(skip(s):end) a(1:skip(s)-1)];   %Shift the amplitude,
   %and compute the mean length.
   surrogate_m(s) = abs(mean(surrogate_amplitude.*exp(1i*p)));
  end
 
  %Compare true mean to surrogates.
  [surrogate_r_mean, surrogate_r_std] = normfit(surrogate_m);
  m_norm_length = (abs(m_raw) - surrogate_r_mean)/surrogate_r_std;  %Create a z-score.
  
  p_values = 1 - normcdf(abs(m_raw), surrogate_r_mean, surrogate_r_std); % Tolga added this part!
