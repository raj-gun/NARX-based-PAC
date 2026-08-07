% This script gives an example of how we computed the performance of the
% PAC methods. It uses the functions that we included in the Supplementary Material
% It is only for "MI" and "Direct Estimate" (modulationindex_directestimate.m), 
% but can be easily changed to the other methods provided as functions GLM (general_linear_index.m) and 
% MI with statistics (modulationindex_stats.m)

% Tolga Ozkurt, tolgaozkurt@gmail.com

clear;clc;

flow1 = 2:2:40;
fhigh1 = 40:10:120;

fs = 1000; % sampling frequency
lowbw = 2;
highbw = 10;

flow2 = flow1+lowbw;
fhigh2 = fhigh1+highbw;

flow = (flow1 + flow2) / 2;
fhigh = (fhigh1 + fhigh2) / 2;

% Compute performance indices 

% range lowfreqs: 12-18 Hz; range highfreqs: 50-100 Hz
indL1 = 6; indL2=9; indH1 = 2; indH2 = 5;
matrlen1 = 25; matrlen2 = 31;

const_coef = (matrlen1 * matrlen2) / ( (indL2 - indL1 + 1) * (indH2 - indH1 + 1) );



%% run 100 iterations

%% find the detection ratio!
% noise levels corresponding to SNR of -10, -5, 0, 10 , 20 dB
noise_levels = [3 1.7 0.93 0.28 0.1];


for noise_ind = 1:5

    clc
    detection_vec1 = zeros(1,100); detection_vec2 = zeros(1,100); 
    
    % average PAC portraits
    DE_vec{noise_ind} = zeros(length(flow),length(fhigh));
    MI_vec{noise_ind} = zeros(length(flow),length(fhigh));

    % 100 iterations with random additive Gaussian noise contaminated
    % simulated data
    for ii = 1:100

        ii
        pacdat = synthesize_pac(noise_levels(noise_ind));
        [DE MI] = modulationindex_directestimate(pacdat, fs,flow1,fhigh1,lowbw,highbw); close all;
        DE_vec{noise_ind} = DE_vec{noise_ind} + DE;
        MI_vec{noise_ind} = MI_vec{noise_ind} + MI;
        
        [a b] = find(DE == max(DE(:))); 
        if (a >=indL1) & (a <=indL2) & (b >= indH1) & (b <= indH2)  % criterion for identification of PAC
            detection_vec1(ii) = 1;
        end


        [a b] = find(MI == max(MI(:))); 
        if (a >=indL1) & (a <=indL2) & (b >= indH1) & (b <= indH2)  
            detection_vec2(ii) = 1;
        end



    end


% this denotes how many times PAC could NOT be detected for DE
nondetection_mat1(noise_ind) = length(find(detection_vec1 == 0)); 
% this denotes how many times PAC could NOT be detected for raw MI 
nondetection_mat2(noise_ind) = length(find(detection_vec2 == 0));

 DE_vec{noise_ind} = DE_vec{noise_ind} ./ 100; % take the average
 MI_vec{noise_ind} = MI_vec{noise_ind} ./ 100;


end

%% Plot the resultant PAC portraits

% direct estimate

figure
subplot(1,5,1)
imagesc(flow,fhigh,DE_vec{1}')
set(gca,'XTick',0:10:40)
axis xy

subplot(1,5,2)
imagesc(flow,fhigh,DE_vec{2}')
set(gca,'XTick',0:10:40)
axis xy

subplot(1,5,3)
imagesc(flow,fhigh,DE_vec{3}')
set(gca,'XTick',0:10:40)
axis xy

subplot(1,5,4)
imagesc(flow,fhigh,DE_vec{4}')
set(gca,'XTick',0:10:40)
axis xy

subplot(1,5,5)
imagesc(flow,fhigh,DE_vec{5}')
set(gca,'XTick',0:10:40)
axis xy

% modulation index

figure
subplot(1,5,1)
imagesc(flow,fhigh,MI_vec{1}')
set(gca,'XTick',0:10:40)
axis xy

subplot(1,5,2)
imagesc(flow,fhigh,MI_vec{2}')
set(gca,'XTick',0:10:40)
axis xy

subplot(1,5,3)
imagesc(flow,fhigh,MI_vec{3}')
set(gca,'XTick',0:10:40)
axis xy

subplot(1,5,4)
imagesc(flow,fhigh,MI_vec{4}')
set(gca,'XTick',0:10:40)
axis xy

subplot(1,5,5)
imagesc(flow,fhigh,MI_vec{5}')
set(gca,'XTick',0:10:40)
axis xy