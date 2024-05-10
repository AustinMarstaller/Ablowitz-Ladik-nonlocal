%% Create a simple signal with two frequencies
dt = .001;
t = 0:dt:1;
f = sin(2*pi*50*t) + sin(2*pi*120*t);
% Sum of 2 frequencies
f = f + 2.5*randn(size(t));

%% Compute the Fast Fourier Transform FFT
n = length(t);
fhat = fft(f,n);
% Compute the fast Fourier transform
PSD = fhat.*conj(fhat)/n;
% Power spectrum (power per freq)
freq = 1/(dt*n)*(0:n-1);
% Create x-axis of frequencies in Hz
L = 1:floor(n/2);
% Only plot the first half of freqs

%% Use the PSD to filter out noise
indices = PSD>100;
% Find all freqs with large power
PSDclean = PSD.*indices;
% Zero out all others
fhat = indices.*fhat;
% Zero out small Fourier coeffs. in Y
ffilt = ifft(fhat);
% Inverse FFT for filtered time signal

%% Plot the PSD which is unfiltered
figure 
subplot(2,1,1)
plot(freq,PSD, 'LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Frequency, Hz')
ylabel('PSD')

%% Plot against wavenumbers instead
kappa = (2*pi)/(n) * [-n/2:n/2-1];
subplot(2,1,2)
plot(kappa,PSD, 'LineWidth',2.5)
set(gca,'FontSize',16)
set(gca,'XTick',-pi:pi/2:pi)
ax.XTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};

xlabel('Wavenumber')
ylabel('PSD')