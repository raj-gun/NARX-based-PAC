function tm_frq_plt(data, Fs, fftn)
%TM_FRQ_PLT Display a signal in time, frequency, and time-frequency domains.
%   The layout combines a time trace, a one-sided magnitude-spectrum view,
%   and a Morse-wavelet scalogram. Linked axes facilitate comparison of
%   temporal events with their spectral content.
%   Inputs
%   ------
%   data  : Signal vector to visualise.
%   Fs    : Sampling frequency in Hz.
%   fftn  : FFT length used for the magnitude-spectrum panel.
%   Output
%   ------
%   This function returns no variables; it creates and formats a figure.

%% Create the time trace, magnitude spectrum, and wavelet scalogram
Ts=1/Fs;
figure;
axis tight
ax1=subplot(2,2,2); plot(data); axis([1,inf,-inf,inf]);
set(gca,'XTickLabel',[]);
ax2=subplot(2,2,3);
w = 0:Fs/fftn:Fs-(Fs/fftn);
plot(abs(Ts.*fft(data,fftn)), w, 'LineWidth', 2); axis([0,inf,1,Fs/2]);
set(gca,'XDir','reverse'); set(gca,"yscale","log"); set(gca, 'FontSize', 18);
set(gca,'XTickLabel',[]); %set(gca,'LineWidth',0.5,'TickLength',[0.05 0.05]);
ax3=subplot(2,2,4); 
[cwt_mat, freqs] = cwt(data, 'morse', Fs);
imagesc(1:length(data), freqs, abs(cwt_mat) ); set(gca,"yscale","log"); shading interp;
set(gca,'LineWidth',1.5,'TickLength',[0.05 0.05]); set(gca, 'FontSize', 18); set(gca,'YDir','normal'); axis([1,inf, 1,Fs/2, 0,inf]);
set(gca,'YTickLabel',[]); set(gca,"yscale","log");
%% Arrange the three axes and link their corresponding dimensions
ax3.Position(3) = ax3.Position(3)*2;
ax3.Position(4) = ax3.Position(4)*2;
ax3.Position(1) = 0.275;
ax3.Position(2) = ax3.Position(2);
ax2.Position(4) = ax3.Position(4);
ax1.Position(3) = ax3.Position(3);
ax2.Position(3) = 0.195;
ax1.Position(4) = 0.2;
ax2.Position(1) = ax3.Position(1)-0.2;
ax1.Position(2) = ax2.Position(2)+0.69;
ax1.Position(1) = ax3.Position(1);
ax2.Position(2) = ax3.Position(2);
linkaxes([ax1,ax3],'x'); linkaxes([ax2,ax3],'y');
end
