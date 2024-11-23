function [power,freq] = powerSpectrum(Signal,Fs)

N = length(Signal);
xdft = fft(Signal);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(Signal):Fs/2;

power = 10*log10(psdx);
end

