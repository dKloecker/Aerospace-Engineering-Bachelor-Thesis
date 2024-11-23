function [fftFrequencies,P1] = signalFFT(Fs,Signal)

L = length(Signal);
fftFrequencies = Fs*(0:(L/2))/L;

Y = fft(Signal);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

end

