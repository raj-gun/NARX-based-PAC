%   adapted from the code written for Kramer et al. (2008) 
%    by Tolga Ozkurt (2011)
%   to implement a General Linear Model (GLM) phase-amplitude coupling estimator

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
%    mod2d_raw1  = The two-dim GLM estimate
%    mod2d_raw2 =  The two-dim GLM estimate from regreession coefffients
%    mod2d_raw3 =  The two-dim spurious term removed robust GLM estimate 
%
%    flow   = The frequency axis for the lower phase frequencies.
%    fhigh  = The frequency axis for the higher amplitude frequencies.
%
%
%

 
function [mod2d_raw1, mod2d_raw2,mod2d_raw3, flow, fhigh] = ...
    general_linear_index(d, srate,flow1,fhigh1,bwlow,bwhigh)
 
  
  flow2 = flow1+bwlow;                                %2-41 Hz with 1 Hz steps.
  
  fhigh2 = fhigh1+bwhigh;                            %5-205 Hz with 5 Hz steps.
  

  mod2d_raw1 = zeros(length(flow1), length(fhigh1));
  mod2d_raw2 = zeros(length(flow1), length(fhigh1));
  mod2d_raw3 = zeros(length(flow1), length(fhigh1));
 
  
  for i=1:length(flow1)
      theta=eegfilt(d,srate,flow1(i),flow2(i));  %Compute the low freq signal.
      theta=theta(srate:length(d)-srate-1);           %Drop the first and last second.
      phase = angle(hilbert(theta));                  %Compute the low freq phase.
      % ['Loops remaining = ', num2str(length(flow1)-i)]
      parfor j=1:length(fhigh1)
        gamma=eegfilt(d,srate,fhigh1(j),fhigh2(j));%Compute the high freq signal.
        gamma=gamma(srate:length(d)-srate-1);         %Drop the first and last second.
        amp = abs(hilbert(gamma));                 %Compute the high freq amplitude.

        %Compute the GLM estimate
        [mi1 mi2 mi3] = general_linear_regress_index(amp, phase);
        mod2d_raw1(i,j) = abs(mi1);
        mod2d_raw2(i,j) = abs(mi2);
        mod2d_raw3(i,j) = abs(mi3);
      end
  end
 
  %Plot the two-dimensional index.
  
  flow = (flow1 + flow2) / 2;
  fhigh = (fhigh1 + fhigh2) / 2;
  
  figure
  imagesc(flow, fhigh, mod2d_raw1');  colorbar;
  axis xy
  set(gca, 'FontSize', 18);
  % xlabel('Phase Frequency (Hz)');  ylabel('Envelope Frequency (Hz)');
  

  figure
  imagesc(flow, fhigh, mod2d_raw2');  colorbar;
  axis xy
  set(gca, 'FontSize', 18);
  % xlabel('Phase Frequency (Hz)');  ylabel('Envelope Frequency (Hz)');
  
  
  figure
  imagesc(flow, fhigh, mod2d_raw3');  colorbar;
  axis xy
  set(gca, 'FontSize', 18);
  % xlabel('Phase Frequency (Hz)');  ylabel('Envelope Frequency (Hz)');
  
 %% end of the function
 
 function [m_raw1, m_raw2, m_raw3] = general_linear_regress_index(a, p)

   % General Linear regression Model PAC estimator
   % a: amplitude, p: phase


numpoints = length(a);

X = [cos(p)' sin(p)' ones(numpoints,1)];
[beta_coef,redundant_term, error_trms] = regress(a', X);
  
% standard GLM (Penny et al., 2008)
m_raw1 = sqrt( (sum(a.*a) - sum(error_trms.*error_trms)) / sum(a.*a) );
 
% equivalent to standard GLM shown by (Ozkurt and Schnitzler, 2011)
m_raw2 = sqrt ( numpoints * (beta_coef(1)^2 + beta_coef(2)^2 + beta_coef(3)^2) / sum(a.*a)  );
 
% robust GLM without the spurious term (beta_3) derived by (Ozkurt and Schnitzler, 2011)  
m_raw3 = 0.5 * sqrt ( (beta_coef(1)^2 + beta_coef(2)^2) / sum(a.*a)  );
 