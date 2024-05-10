uhat = fft(uu, L);
PSD = uhat.*conj(uhat)/L;

kappa = (2*pi/L)*[-L/2:L/2-1]; % Discrete wavenumbers?

plot(kappa, PSD(end,:))